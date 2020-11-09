###################################
#
# utility script
#
# written by M. Fidino 1/6/2020 mdy
#
###################################

###########
# update_foia_columns
###########

# Some general quality of life naming improvements
update_foia_columns <- function(data){
	colnames(data) <- c(
		"type",
		"request_number",
		"date",
		"year",
		"month",
		"block",
		"description")
	return(data)
}
##################
# drop_uncertain_species
##################

# This will run you through the records that have some uncertainty
#   associated to what species the record is actually about.
#   You then have to read through each and enter input on wheter
#   to keep it.
drop_uncertain_species <- function(data, species){
	
	# create a vector to determine which records to trust
	#  All start out untrusted
	to_keep <- rep(
		FALSE,
		nrow(data)
	)
	# These are regex bits that indicate
	#  that a record should be reviewed
	uncertainties <- c(
		paste0(species, "s?\\s+or"),
		paste0("\\sor\\s+",species),
		paste0("\\sor\\s+a\\s+",species),
		"unsure",
		"uncertain",
		paste0(species,"\\s?/"),
		paste0("/\\s?",species),
		"unknown",
		"don't know",
		"do not know",
		"dont know"
	)
	
	# combine them
	uncertain_regex <- paste(
		uncertainties,
		collapse = "|"
	)
	
	# flag the issues
	flagged <- grep(
		uncertain_regex,
		data$description
	)
	
	to_keep[-flagged] <- TRUE
	
	# loop through and determine which ones to keep
	cat(paste(length(flagged), "records to review\n"))
	for(i in 1:length(flagged)){
		tmp <- utils::menu(
			choices = c(TRUE,	FALSE),
			title = data$description[flagged[i]]
		)
		# if the first option, keep it.
		to_keep[flagged[i]] <- ifelse(
			tmp == 1,
			TRUE,
			FALSE
		)
	}
	# keep the correct records
	data <- data[to_keep,]
	
	# write the data
	write.csv(data, paste0("./data/conflict_clean/", species,".csv"),
						row.names = FALSE)
} 

##########################
# geocode_wrapper
##########################

# A convience function to geocode data from the google api
#   This is just so it can automatically run the queries each
#   day at 6 AM.
geocode_wrapper <- function(data, address_column){
	# the number of records to geocode
	ndata <- nrow(data)
	# figure out number of days this will take at 2.5 K per day
	ndays <- ceiling(ndata / 2500)
	# spit out report
	cat(paste("This will take", ndays, "days to geocode these data...\n"))
	cat("Starting geocode process...\n")
	# split into batches
	code_groups <- factor((1:nrow(data) %/% 2500))
	data <- split(data, f = code_groups)
	batch <- 1
	while(batch < ndays){
		# do the geocoding
		ans <- ggmap::geocode(as.character(data[[batch]][,address_column]))
		data$lon <- ans$lon
		data$lat <- ans$lat
		# get current time
		time <- lubridate::ymd_hms(Sys.time())
		# get 6 AM tomorrow
		tomorrow <- paste(
			lubridate::today() + 1,
			"06:00:00 UTC"
		)
		# calculate time to the next day
		time_to_tomorrow <- difftime( 
			lubridate::ymd_hms(
				tomorrow
			),
			time,
			units = "secs"
		)
		# Go to sleep if more to do, stop this while statement by
		#  incrementing batch to > ndays.
		if(batch < ndays){
			cat(paste("Batch", batch,"of",ndays, "complete.",
								"Going to sleep until", tomorrow,"..."))
			Sys.sleep(as.numeric(time_to_tomorrow))
			batch <- batch + 1
		} else {
			cat(paste("Batch", batch,"of",ndays, "complete."))
			batch <- batch + 1
		}
	}
	# return the data
	to_return <- dplyr::bind_rows(data)
	return(to_return)
}

######################
# point_count
######################

# This function takes in a raster and a set of coordinates
point_count <- function(my_raster, my_pts){
	if(!inherits(my_raster, "Raster")){
		stop("my_raster must be a Raster* object.")
	}
	if(!"sf" %in% class(my_pts)){
		stop("my_pts must be a sf object.")
	}
	# Collect the spatial coordinates
	my_pts <- sf::st_coordinates(my_pts)
	# make a raster of zeroes like the input
	tmp <- my_raster
	tmp[] <- 0
	# get the cell index for each point and make a table:
	counts <- table(
		raster::cellFromXY(tmp, my_pts)
	)
	# assign the counts to the raster
	tmp[as.numeric(names(counts))] <- counts
	# collect the values from the raster
	
	return(tmp)
}

#########################
# order seasons
#########################


# Seasons are ordered like FA13, SP13, etc. They need to get ordered
#  chronologically, not alphabetically
order_seasons <- function(x){
	unq_seasons <- unique(x)
	month <- substr(unq_seasons,1,2)
	
	tmp <- sapply(month, function(y) switch(y, "WI" = 0.1,
																					"SP" = 0.2,
																					"SU" = 0.3,
																					"FA" = 0.4)
	)
	
	year <- as.numeric(paste0("20", substr(unq_seasons, 3,4)))
	
	to_order <- year + tmp
	
	unq_seasons <- unq_seasons[order(to_order)]
	
	return(unq_seasons)
	
}

# There are a few sites that should be within the city limits but are not
#  This function pushes them to the nearest cell. A lot of this was 
#  borrowed from rSDM::points2nearestcell(), though that package
#  is not currently available in R version 4.0.3. 

# locs = sf data.frame of spatial points
# r = a single raster layer
points_to_nearest_raster_cell <- function(
	locs , r, distance = NULL ) {
	pid <- raster::cellFromXY(
		r,
		sf::st_coordinates(locs)
	)
	
	miss <- which(is.na(values(r)[pid]))
	
	## if there are NA cells...
	
	if (sum(miss) > 0) {
		# points without data
		coord_miss <- sf::st_coordinates(locs[miss, ])  
		# get coordinates of cells with data
		if (nlayers(r) > 1){
			r <- raster::raster(r, 1)
		}
		cells_notNA <- raster::rasterToPoints(r, spatial = TRUE)  
		coord_ras <- sp::coordinates(cells_notNA)
		cell_id <- factor(seq_len(nrow(coord_ras)))
		
		# find the nearest raster cell for each point with missing data
		nearest_cell <- class::knn1(coord_ras, coord_miss, cell_id)
		
		new_coords <- matrix(coord_ras[nearest_cell, ], ncol = 2)
		colnames(new_coords) <- c("X", "Y")
		
		if (!is.null(distance)){
			
			# calculate distances between old and new coordinates
			distances <- raster::pointDistance(
				coord_miss, 
				new_coords,
				lonlat = raster::isLonLat(locs)
			)
			# if distance below threshold, accept, otherwise keep old coordinates
			x <- ifelse(distances < distance, new_coords[,1], coord_miss[,1])
			y <- ifelse(distances < distance, new_coords[,2], coord_miss[,2])
			new_coords <- cbind(X = x, Y = y)
			
		}
		
		
		# assign new coordinates to those points
		og_coords <- sf::st_coordinates(locs)
		og_coords[miss, ] <- new_coords
		
		new_geometry <- sf::st_as_sf(
			data.frame(og_coords),
			coords = c(1,2)
		)
		sf::st_geometry(locs) <- sf::st_geometry(new_geometry)
		
	} else message("All points fall within a raster cell")
	
	return(locs)
	
}
