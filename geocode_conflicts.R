####################################
#
# Geocoding human-wildlife conflicts
#
# Written by M. Fidino 1/6/2020  mdy
#
####################################

# Note:
#   In the event that you are querying over 2.5K locations this
#   script will run over multiple days so you do not incur
#   charges from google. As such, once you run this it's
#   best to just leave it in the background. 

# load packages
library(ggmap)
library(keyring)
library(dplyr)

# load utility script to bring some helper functions
source("utility_script.R")

# I've stored my api key in the keyring package so I don't
#  store it in a text file. You can uncomment this following
#  code once you have your google api set up

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

# Raccoon

# bring in the clean data
raccoon <- read.csv(
	"./data/conflict_clean/raccoon.csv",
	stringsAsFactors = FALSE
)

# Add Chicago, Illinois to the block the conflict occurred.
raccoon$block <- paste0(
	raccoon$block,
	", Chicago, Illinois"
)

# Use wrapper function I wrote for ggmap::geocode
raccoon_geocoded <- geocode_wrapper(
	data = raccoon,
	address_column = "block"
)

# save as a csv
write.csv(
	raccoon_geocoded,
	"./data/conflict_clean/raccoon.csv",
	row.names = FALSE
)
	


