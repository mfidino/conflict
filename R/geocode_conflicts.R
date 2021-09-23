####################################
#
# Geocoding human-wildlife conflicts
#
# Written by M. Fidino 1/6/2020  mdy
#
####################################


# load packages
library(ggmap)
library(keyring)
library(dplyr)

# load utility script to bring some helper functions
source("./R/sourcer.R")

species <- "coyote"

# I've stored my api key in the keyring package so I don't
#  store it in a text file. You can uncomment this following
#  code once you have your google api set up
# Assuming you have already set up the API, I went here
# https://console.cloud.google.com/google/maps-apis/credentials

########################################
# UNCOMMENT AND FILL IN YOUR API BELOW
#  TO SAVE YOUR GOOGLE API.
#
## Create your keyring
# keyring::keyring_create("hw_conflict")
#
## Attach your key to the keyring
# keyring::key_set_with_value(
# 	service = "maps_api",
# 	password = "YOUR GOOGLE API",
# 	keyring = "hw_conflict"
# )
########################################

## With the above code you can then call your API with
# keyring::key_get(
#   service = "maps_api",
#   keyring = "hw_conflict"
# )

# register with google
ggmap::register_google(
	key = keyring::key_get(
		service = "maps_api",
		keyring = "hw_conflict"
	),
	account_type = "standard"
)


# bring in the clean data
my_data <- read.csv(
	paste0("./data/conflict_clean/",species,".csv"),
	stringsAsFactors = FALSE
)

# Add Chicago, Illinois to the block the conflict occurred.
my_data$block <- paste0(
	my_data$block,
	", Chicago, Illinois"
)

# Use wrapper function I wrote for ggmap::geocode
my_data_geocoded <- geocode_wrapper(
	data = my_data,
	address_column = "block"
)


# Check to see if there are any NA coordinates
if(any(is.na(my_data_geocoded$lon))){
  tmp <- geocode_wrapper(
  	my_data[is.na(my_data_geocoded$lon),],
  	address_column = "block"
  )
  tmp <- data.frame(tmp)
  my_data_geocoded[is.na(my_data_geocoded$lon),c("lon", "lat")] <- 
  	tmp[,c("lon", "lat")]
}

# save as a csv
write.csv(
	my_data_geocoded,
	paste0("./data/conflict_clean/", species,".csv"),
	row.names = FALSE
)
	

