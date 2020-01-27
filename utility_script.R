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
