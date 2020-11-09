##########################
#
# Formatting data to be fit by model
#
# Written by M. Fidino 1/16/2020
#
#########################

species <- "raccoon"
source("sourcer.R")


packs <- c("lubridate", "raster", "sp", "sf")

package_load(packs)

# Set up the presence only data

# Step 1. Split into different time periods.
podat <- read.csv(
	paste0("./data/conflict_clean/", species,".csv"),
	stringsAsFactors = FALSE
)

# create date object, but remove any records that lack date info
if(any(is.na(podat$date) | podat$date == "")){
	podat <- podat[-which(is.na(podat$date) | podat$date == ""),]
}
podat$date <- lubridate::mdy(podat$date)
podat$year <- lubridate::year(podat$date)
podat$month <- lubridate::month(podat$month)

# Group podat into four different seasons (and by year)
#   1. Dec, Jan, and Feb
#   2. Mar, Apr, May
#   3. Jun, Jul, Aug
#   4. Sep, Oct, Nov

# get range of years
uyear <- unique(podat$year)

# construct groups based on years, we need the year preceeding
#  to more easily construct the initial year of sampling
my_years <- c(min(uyear) - 1, uyear)

# This will make a set of groups based on the year and month,
#  removing the first 11 (Jan through December) as we don't
#  need those for the 'first' year.
year_month <- expand.grid(
	list(
		month = 1:12,
		year = my_years
	)
)[-c(1:11),]

# Drop the 12 2013 row because we don't have any pa data for that
year_month <- year_month[-nrow(year_month),]

# calculate the number of seasons
n_season <- ceiling(nrow(year_month) / 3)
# and from there construct the season grouping. This will
#   create 3 values for each season, but then the last season
#   may not have 3 months available. So we subset this down
#   to the number of rows in year_month.
year_month$season <- rep(seq(1:n_season), each = 3)[1:nrow(year_month)]

# combine that with presence only data
podat <- dplyr::inner_join(
	podat,
	year_month,
	by = c("year", "month")
)

# sort by season
podat <- podat[order(podat$season),]

# The next thing we need to do is get a count of conflicts across our raster.
chicago_raster <- readRDS("./data/all_raw_layers.rds")

# make podat spatial
podat <- st_as_sf(
	podat,
	coords = c("lon", "lat"),
	crs = 4326
)
# then convert to UTM (the Chicago raster)
podat <- st_transform(
	podat,
	crs = crs(chicago_raster)
)

# We need to get the count in cells, point_count()
#  below does this for us (sourced at top of script).
#  We also need to do it for each
#  season
my_counts <- vector(
	"list",
	length = n_season
)
for(i in 1:n_season){
	if(sum(podat$season == i) == 0){
		next
	}
	tmp_podat <- podat[which(podat$season == i),]
  my_counts[[i]] <- point_count(
  	chicago_raster$imperv,
  	tmp_podat
  )
  names(my_counts[[i]]) <- species
  
  # Remove any counts outside of Chicago
  #  and also turn 0 to NA
  tmp <- values(my_counts[[i]][[species]])
  tmp[tmp == 0] <- NA
  tmp[which(is.na(values(chicago_raster$imperv)) & tmp > 0)] <- NA
  values(my_counts[[i]]) <- tmp
  rm(tmp)
  rm(tmp_podat)
}

id_vals <- 1:sum(!is.na(values(chicago_raster$canopy)))
chicago_raster$id <- NA
id_to_fill <- values(chicago_raster$id)
id_to_fill[!is.na(values(chicago_raster$canopy))] <- id_vals
values(chicago_raster$id) <- id_to_fill


# We now have a PO raster for each season of sampling.
#  we'll need to convert it to a ragged array.

# Create the list object we are going to store our
#  counts in. The Presence only (po) pixel_id
po_pixel_id <- vector(
	"list",
	length = n_season
)
# Loop through each season.
for(i in 1:n_season){
	# Collect the values
  tmp <- values(my_counts[[i]][[species]])
  tmp_ids <- values(chicago_raster$id)

  # determine where they occur in the giant raster
  po_location <- which(tmp > 0)

  # We have multiple records in some squares, so these need to be 
  #  duplicated for the model.
  times_to_rep <- tmp[po_location]
  
  # determine their ID based on only the cells that have values
  po_id <- tmp_ids[po_location]
  
  # Put the presence only location into our list object.
  po_pixel_id[[i]] <- rep(
	  po_id,
	  times = times_to_rep
  )
}

# get the npo count per season
npo_count <- sapply(po_pixel_id, length)

rm(times_to_rep)
rm(tmp)
rm(my_years)
rm(uyear)
rm(podat)

# Step 2. Now we need to pull in the appropriate species data
#  from the detection / non-detection data.

# Bring in and combine all the data

padat <- read.csv(
	"./data/camera_trap_detections_sp10_sp13.csv"
)

# drop 2010 data
padat <- padat[-grep("10", padat$Season),]

padat <- padat %>% 
	dplyr::group_by(
	  StationID, Season
  ) %>% 
	dplyr::summarise(
		J = sum(!is.na(Coyote)),
		Coyote = sum(Coyote, na.rm = TRUE),
		Opossum = sum(Opossum, na.rm = TRUE),
		Raccoon = sum(Raccoon, na.rm = TRUE),
		Redfox = sum(Redfox, na.rm = TRUE),
		Skunk = sum(Skunk, na.rm = TRUE)
	)


padat <- tidyr::pivot_longer(
	padat,
	c("Coyote", "Opossum", "Raccoon", "Redfox", "Skunk"),
	"Species",
	values_to = "count"
)

#  Summer 2013 data

padat2 <- read.csv(
	"./data/summer_2013.csv"
)

padat2$Season <- "SU13"

colnames(padat2) <- c("Species", "StationID", "count", "J", "Season")

padat3 <- read.csv(
	"./data/fall_13.csv"
)
padat3$Season <- "FA13"

colnames(padat3) <- colnames(padat2)

padat <- dplyr::bind_rows(
	list(padat, padat2, padat3)
)

padat$Season <- factor(
	padat$Season,
	order_seasons(padat$Season)
)

# order data in a logical way
padat <- padat[
	order(
		padat$Season,
		padat$Species,
		padat$StationID
	),
]


station_coords <- read.csv(
	"./data/station_coords.csv"
)

sc_spatial <- sf::st_as_sf(
station_coords,
coords = c("Easting", "Northing"),
crs = 32616
)


sc_spatial <- sf::st_transform(
	sc_spatial,
	sf::st_crs(chicago_raster)
)

city_outline <- sf::st_read(
	"D:/GIS/IL City Shape Files", 
	layer = "tl_2017_17_place"
) %>% 
	dplyr::select("NAME")

# get Chicago
chicago <- city_outline[city_outline$NAME == "Chicago",]

# Transform to same projection as the raster
chicago <- sf::st_transform(
	chicago,
	crs = sf::st_crs(
		chicago_raster
	)
)

# crop the Chicagoland region to just Chicago

# The crop only get's us to the rectangular extent of the chicago
#   polygon. We want to NA the values outside of the Chicago polygon
#   Adding a buffer of 1000m to account for sites or conflicts
#   on the edge.

sc_spatial <- sf::st_intersection(
	sc_spatial,
	chicago
)

# There are 5 sites that should be within the city (along the edge of 
#  the raster). We are going to move them to the nearest cell with
#  information.

sc_spatial <- points_to_nearest_raster_cell(
	sc_spatial,
	chicago_raster
)

# add a numeric ID to each site, first we need to get the cell values
#   from the whole raster
pa_full_pixelID <- raster::cellFromXY(
	chicago_raster$imperv,
	sf::st_coordinates(sc_spatial)
)

# And then assign them the actual ID, which is based on the raster cells
#  with information
sc_spatial$pixelID <- values(chicago_raster$id)[pa_full_pixelID]


# convert back to a data.frame
sc_spatial <- data.frame(
	StationID = sc_spatial$StationID,
	pixelID = sc_spatial$pixelID
)

station_coords <- dplyr::inner_join(
	station_coords,
	sc_spatial,
	by = "StationID"
)

# join camera trap data to station coordinates
padat <- dplyr::inner_join(
	padat,
	station_coords,
	by = "StationID"
)

padat <- sf::st_as_sf(
	padat,
	coords = c("Easting", "Northing"),
	crs = 32616
)

padat <- sf::st_transform(
	padat,
	crs = crs(chicago_raster)
)


# filter down to the target species for this analysis
padat$Species <- tolower(padat$Species)
padat <- padat[padat$Species == species,]

# See if any sites need to be dropped because they have no data, this
#  occured because the data we have here is a bit older than all the
#  station coordinates we ahve
drop_these_sites <- padat %>% 
	dplyr::group_by(StationID) %>% 
	dplyr::summarise(dsamp = sum(J)) %>% 
	dplyr::filter(dsamp == 0)

padat <- padat[-which(padat$StationID %in% drop_these_sites$StationID),]

# the ymatrix needs to be in wide format, first column is pixels
y_widen <- data.frame(
	pixelID = padat$pixelID,
	Season = padat$Season,
	count = padat$count
)

y_widen <- tidyr::pivot_wider(
	y_widen,
	names_from = Season,
	values_from = count
)

# same thing with the J matrix
j_widen <- data.frame(
	pixelID = padat$pixelID,
	Season = padat$Season,
	J = padat$J
)
which(is.na(j_widen), arr.ind = TRUE)

j_widen <- tidyr::pivot_wider(
	j_widen,
	names_from = Season,
	values_from = J
)

# a few of these are NA, which means no sampling
j_widen[is.na(j_widen)] <- 0

# go thorugh and put an NA if we did not sample
to_NA <- which(
	j_widen == 0,
	arr.ind = TRUE
)

for(i in 1:nrow(to_NA)){
	tmp <- to_NA[i,]
	# stop if we are overwriting actual data for some reason
	if(y_widen[tmp[1], tmp[2]] > 0 & !is.na(y_widen[tmp[1], tmp[2]])){
		stop()
	}
	y_widen[tmp[1], tmp[2]] <- NA
}

# Get the number of unique cells (i.e., not NA cells in our raster)

G <- sum(
	!is.na(
		raster::values(chicago_raster$canopy)
	)
)

# Scale the raster data
chicago_scaled <- raster::scale(
	chicago_raster
)

# check for correlation among covariates
to_corr <- values(chicago_scaled)
to_corr <- to_corr[complete.cases(to_corr),]

corred <- round(
	cor(to_corr),
	2
)

# get the occupancy covariates
occ_covs <- cbind(1, values(chicago_scaled))

# and drop out all NA values
occ_covs <- occ_covs[complete.cases(occ_covs),
										 -which(colnames(occ_covs) == "id")]

# looks like only grass cover and impervious covar are strongly and
#  negatively correlated (-0.62). The next would be canopy cover and
#  impervious cover (-0.43). Everything else is quite low in terms
#  of correlation. As such, I'm just going to let the Bayesian
#  LASSO take care of this.

my_data <- list(
	# Total number of grid points
	G = G, 
	# The occupancy covariates of the latent state Poisson process
	occ_covs = occ_covs,
	# The detection covariates for the presence only data
	po_det_covs = occ_covs,
	# The detection covariates for the presence / absence data
	pa_det_covs = occ_covs,
	# The pixels that the presence absence data occurred
	pa_pixel = y_widen$pixelID,
	# The pixels that the presence only data occurred							
	po_pixel = unlist(po_pixel_id),
	# What season each opportunistic data point came from
	opp_year = rep(1:length(po_pixel_id), times = npo_count),
	# The number of presence only data points per season
	npo = as.numeric(npo_count),
	# The total number of all presence only data points
	all_npo = sum(npo_count),
	# The number of days a species was detected per site / season
	#   for the presence absence data
	y_pa = as.matrix(y_widen[,-1]),
	# The number of sites presence absence data was sampled
	npa = nrow(y_widen),
	# The log cell area, logged as the parameters are on the log scale
	#   in the model.
	cell_area = log(prod(res(chicago_raster))),
	# Ones for the 'ones' trick in JAGS (for coding up likelihood)
	ones = rep(1, sum(npo_count)),
	# A big constant value for 'ones' trick.
	CONSTANT = 10000,
	# Number of latent parameters
	nlatent = ncol(occ_covs),
	# Number of observational parameters, presence only
	nobs_po = ncol(occ_covs),
	# Number of observational parameters, presence / absence
	nobs_pa = ncol(occ_covs),
	# Number of seasons sampled
	nyear = length(po_pixel_id),
	J = as.matrix(j_widen[,-1])
)


inits <- function(chain){
	gen_list <- function(chain = chain){
		list( 
			z = matrix(1, ncol = my_data$nyear, nrow = my_data$G),
			beta_occ_fill = rnorm(my_data$nlatent-1),
			beta_observation = rnorm(my_data$nobs_pa-1),
			beta_po_fill = rnorm(my_data$nobs_po-1),
			psi_mu = rnorm(1),
			pa_mu = rnorm(1),
			po_mu = rnorm(1),
			lambda_beta_occ = runif(1, 1, 3),
			lambda_pa_det = runif(1, 1, 3),
			lambda_po_det = runif(1, 1, 3),
			psi_tau_season = rgamma(1, 1, 1),
			pa_tau_season = rgamma(1, 1, 1),
			po_tau_season = rgamma(1, 1, 1),
			psi_season = rnorm(my_data$nyear),
			pa_season = rnorm(my_data$nyear),
			po_season = rnorm(my_data$nyear),
			.RNG.name = switch(chain,
												 "1" = "base::Wichmann-Hill",
												 "2" = "base::Marsaglia-Multicarry",
												 "3" = "base::Super-Duper",
												 "4" = "base::Mersenne-Twister",
												 "5" = "base::Wichmann-Hill",
												 "6" = "base::Marsaglia-Multicarry",
												 "7" = "base::Super-Duper",
												 "8" = "base::Mersenne-Twister"),
			.RNG.seed = sample(1:1e+06, 1)
		)
	}
	return(switch(chain,           
								"1" = gen_list(chain),
								"2" = gen_list(chain),
								"3" = gen_list(chain),
								"4" = gen_list(chain),
								"5" = gen_list(chain),
								"6" = gen_list(chain),
								"7" = gen_list(chain),
								"8" = gen_list(chain)
	)
	)
}

m1 <- run.jags(model = "integrated_pp_dynamic.R", 
							 data = my_data, 
							 n.chains = 4, 
							 inits = inits, 
							 monitor = c(
							 	"beta_occ_fill", "beta_observation", "beta_po_fill",
							 	"lambda_beta_occ", "lambda_pa_det", "lambda_po_det",
							 	"psi_mu", "pa_mu", "po_mu", "psi_season",
							 	"pa_season", "po_season", "psi_sd_season",
							 	"pa_sd_season", "po_sd_season"
							 	), 
							 adapt = 1000, 
							 burnin = 2000, 
							 sample = 3000,
							 thin = 5,
							 method = 'parallel',
							 summarise = FALSE)
