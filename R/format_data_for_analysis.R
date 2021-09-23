##########################
#
# Formatting data to be fit by model
#
# Written by M. Fidino 1/16/2020
#
#########################


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

# get only the months where we sample
#year_month <- year_month[seq(2, nrow(year_month), by = 3),]

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
tmp_crs <- crs(podat)
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
  tmp <- raster::values(my_counts[[i]][[species]])
  tmp[tmp == 0] <- NA
  tmp[which(is.na(raster::values(chicago_raster$imperv)) & tmp > 0)] <- NA
  raster::values(my_counts[[i]]) <- tmp
  rm(tmp)
  rm(tmp_podat)
}

id_vals <- 1:sum(!is.na(raster::values(chicago_raster$canopy)))
chicago_raster$id <- NA
id_to_fill <- raster::values(chicago_raster$id)
id_to_fill[!is.na(raster::values(chicago_raster$canopy))] <- id_vals
raster::values(chicago_raster$id) <- id_to_fill


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
  tmp <- raster::values(my_counts[[i]][[species]])
  tmp_ids <- raster::values(chicago_raster$id)

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
	  times = 1
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
		Skunk = sum(Skunk, na.rm = TRUE),
		.groups = "drop_last"
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
# change name of opossum in this dataset
padat2$CommonName <- gsub(
	"Virginia opossum",
	"Opossum",
	padat2$CommonName
)
# as well as skunk
padat2$CommonName <- gsub(
	"Striped Skunk",
	"Skunk",
	padat2$CommonName
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
	"./data/IL City Shape Files", 
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
sc_spatial$pixelID <- raster::values(chicago_raster$id)[pa_full_pixelID]


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
	data.frame(padat),
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

# get the coordinates
ccord <- raster::coordinates(chicago_scaled)

ccord <- st_as_sf(
	data.frame(ccord),
	coords = c("x", "y"),
	crs = 32616
)
ccord <- sf::st_transform(
	ccord,
	crs = 4326
)
ccord <- sf::st_coordinates(ccord)
# check for correlation among covariates
to_corr <- raster::values(chicago_scaled)
to_corr <- to_corr[complete.cases(to_corr),]

corred <- round(
	cor(to_corr),
	2
)

# get the occupancy covariates

occ_covs <- raster::values(chicago_scaled)

# and drop out all NA values
ccord <- ccord[complete.cases(occ_covs),]
occ_covs <- occ_covs[complete.cases(occ_covs),
										 -which(colnames(occ_covs) == "id")]




# looks like only grass cover and impervious cover are strongly and
#  negatively correlated (-0.62). The next would be canopy cover and
#  impervious cover (-0.43), and then imperv and houses (0.32).
#  Everything else is quite low in terms
#  of correlation. Regardless, we need to do some amount of dimension
#  reduction.

 var_pca <- prcomp(
 	occ_covs[,c("canopy", "grass", "imperv", "houses")],
 	center = FALSE
 	)
 
 # Here are the loadings from this analysis
 # canopy  0.2784312 -0.8189926  0.2674786 -0.4244789
 # grass   0.5682957  0.2921264 -0.5502036 -0.5375669
 # imperv -0.6672976  0.1781493  0.1044941 -0.7155821
 # houses -0.3927213 -0.4606257 -0.7840996  0.1370465
 
 # multiply everything by -1 so positive generally means more urban
 var_pca$x <- var_pca$x * -1
 var_pca$rotation <-  var_pca$rotation * -1
 
 # looking better
 #          PC1   PC2   PC3   PC4
 # canopy -0.28  0.82 -0.27  0.42
 # grass  -0.57 -0.29  0.55  0.54
 # imperv  0.67 -0.18 -0.10  0.72
 # houses  0.39  0.46  0.78 -0.14
 
 # first two explain 75% of the variability in the data.
 # URB1 = negative is more canopy and grass, positive is more imperv and houses
 # URB2 = negative is more grass and imperv, positive is more houses and canopy
 
 # we'll add those in and check for correlations again
 
 occ_covs <- cbind(var_pca$x[,1:2], occ_covs[,c("income", "vacancy")])
 
 # generate the spatial spline stuff
 tmp_dat <- cbind(1, ccord, occ_covs)
 colnames(tmp_dat)[1:3] <- c("y", "E", "N") 
 tmp_dat <- data.frame(tmp_dat)
 
 
 jags.file <- "./JAGS/test.jags"

 # Thi
 offie <- log(prod(res(chicago_raster)/100))
 gam_dat <- jagam(y ~ s(E,N, k = 10, bs = "ds", m = c(1,0.5)),
 							data = tmp_dat, file = jags.file, 
 							family = "binomial")
 # change initial value of model intercept
 gam_dat$jags.ini$b[1] <- -3.2
 
 psi_covs <- occ_covs[,c("PC1", "income", "vacancy")]
 po_det_covs <-  occ_covs[,c("PC2", "income", "vacancy")]
 
 # reduce detection covaraites for pa data to just pc1 and pc2
 pa_det_covs <- occ_covs[,c("PC1", "PC2")]
 
 my_data <- list(
 	# Total number of grid points
 	G = G, 
 	# The occupancy covariates of the latent state Poisson process
 	occ_covs = psi_covs,
 	# The detection covariates for the presence only data
 	po_det_covs = po_det_covs,
 	# The detection covariates for the presence / absence data
 	pa_det_covs = pa_det_covs,
 	# The pixels that the presence absence data occurred
 	pa_pixel = y_widen$pixelID,
 	# The pixels that the presence only data occurred							
 	po_pixel = unlist(po_pixel_id),
 	# What season each opportunistic data point came from
 	opp_year = rep(1:(length(po_pixel_id)), times = npo_count),
 	# The number of presence only data points per season
 	npo = as.numeric(npo_count),
 	# The total number of all presence only data points
 	all_npo = sum(npo_count),
 	# The number of days a species was detected per site / season
 	#   for the presence absence data
 	y_pa = as.matrix(y_widen[,-c(1)]),
 	# The number of sites presence absence data was sampled
 	npa = nrow(y_widen),
 	# The log cell area, logged as the parameters are on the log scale
 	#   in the model.
 	cell_area = log(prod(res(chicago_raster)/100)),
 	# Ones for the 'ones' trick in JAGS (for coding up likelihood)
 	ones = rep(1, sum(npo_count)),
 	# A big constant value for 'ones' trick.
 	CONSTANT = 10000,
 	# Number of latent parameters
 	nlatent = ncol(psi_covs),
 	# Number of observational parameters, presence only
 	nobs_po = ncol(po_det_covs),
 	# Number of observational parameters, presence / absence
 	nobs_pa = ncol(pa_det_covs),
 	# Number of seasons sampled
 	nyear = length(po_pixel_id),
 	J = as.matrix(j_widen[,-c(1)]),
 	# GAM stuff
 	S1 = gam_dat$jags.data$S1,
 	zero = gam_dat$jags.data$zero,
 	X = gam_dat$jags.data$X,
 	nspline = length(gam_dat$jags.data$zero)
 )

  y_pa_init <- my_data$y_pa
  y_pa_init[is.na(y_pa_init)] <- 0
  y_pa_init[!is.na(my_data$y_pa)] <- NA

 my_inits <- function(chain){
	gen_list <- function(chain = chain){
		list(
			z = matrix(1, ncol = my_data$nyear, nrow = my_data$G),
			beta_occ = rnorm(my_data$nlatent, 0, 0.25),
			b = matrix(
				rnorm(my_data$nspline * my_data$nyear,
							rep(
								gam_dat$jags.ini$b,
								each = my_data$nyear),
							0.1),
				nrow = my_data$nspline,
				ncol = my_data$nyear),
			lambda_gam = rgamma(1, 1,1),
			gam_tau = rgamma(1,1,1),
			beta_pa_det = rnorm(my_data$nobs_pa, 0, 0.25),
			beta_po_det = rnorm(my_data$nobs_po, 0, 0.25),
			beta_po_det_mu = rnorm(my_data$nobs_po),
			beta_po_det_tau = rgamma(my_data$nobs_po,1,1),
			psi_mu = rnorm(1, -5, 0.25),
			pa_mu = rnorm(1, -2.75, 0.25),
			po_mu = rnorm(1, 3, 0.25),
			lambda_beta_occ = runif(1, 1, 2),
			lambda_pa_det = runif(1, 1, 2),
			lambda_po_det = runif(1, 1, 2),
			psi_tau_season = rgamma(1, 2, 4),
			pa_tau_season = rgamma(1, 2, 4),
			po_tau_season = rgamma(1, 2, 4),
			psi_season = rnorm(my_data$nyear, 0, 0.25),
			pa_season = rnorm(my_data$nyear, 0, 0.25),
			po_season = rnorm(my_data$nyear,0, 0.25),
			y_pa = y_pa_init,
			.RNG.name = switch(chain,
												 "1" = "base::Wichmann-Hill",
												 "2" = "base::Wichmann-Hill",
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

 