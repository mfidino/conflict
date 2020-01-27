##########################
#
# Formatting data to be fit by model
#
# Written by M. Fidino 1/16/2020
#
#########################

species <- "raccoon"
source('utility_script.R')

library(lubridate)
library(raster)
library(sp)
library(sf)
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
# calculate the number of seasons
n_season <- ceiling(nrow(year_month) / 3)
# and from there construct the season grouping. This will
#   create 3 values for each season, but then the last season
#   may not have 3 months available. So we subset this down
#   to the number of rows in year_month.
year_month$season <- rep(seq(1:n_season), each = 3)[1:nrow(year_month)]

# combine that with presence only data
podat <- dplyr::left_join(
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

# We need to get the count in cells, the function
#  below does this for us. We also need to do it for each
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

  # determine where they occur
  po_location <- which(tmp > 0)

  # We have multiple records in some squares, so these need to be 
  #  duplicated for the model.
  times_to_rep <- tmp[po_location]
  
  # Put the presence only location into our list object.
  po_pixel_id[[i]] <- rep(
	  po_location,
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



