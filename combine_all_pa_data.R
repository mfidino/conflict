# combining all of the presence / absence data together

library(dplyr)
library(tidyr)

d1 <- read.csv(
	"./data/camera_trap_detections_sp10_sp13.csv",
	stringsAsFactors = FALSE
)


test <- d1 %>% 
	group_by(StationID, SurveyID, Season) %>% 
	summarise(Coyote = sum(Coyote, na.rm = TRUE),
						Opossum = sum(Opossum, na.rm = TRUE),
						Redfox = sum(Redfox, na.rm = TRUE),
						Raccoon = sum(Raccoon, na.rm = TRUE),
						Skunk = sum(Skunk, na.rm = TRUE))

test2 <- d1 %>% 
	group_by(StationID, SurveyID, Season) %>% 
	summarise(J = length(Coyote) - sum(is.na(Coyote)))

test[which(test2$J == 0),4:8] <- NA

test$Season <- factor(test$Season,
											levels = c("SP10", "SU10", "FA10", "WI11",
																 "SP11", "SU11", "FA11", "WI12",
																 "SP12", "SU12", "FA12", "WI13",
																 "SP13", "SU13", "FA13"))

test2$Season <- factor(test2$Season, levels = levels(test$Season))

test <- test[order(test$Season, test$SurveyID),]
test2 <- test2[order(test2$Season, test2$SurveyID),]
test$J <- test2$J

test3 <- pivot_longer(
					test,
					cols = Coyote:Skunk,
					names_to = "CommonName",
					values_to = "count"
					)

test3 <- test3[order(test3$Season, test3$CommonName, test3$StationID),]

test3 <- test3[,c("CommonName", "StationID", "Season", "count", "J")]


# bring in summer 2013
su13 <- read.csv(
	"./data/summer_2013.csv",
	stringsAsFactors = FALSE
)


su13$Season <- "SU13"

su13 <- su13[,c(1,2,5,3,4)]
colnames(su13) <- colnames(test3)

su13$CommonName <- gsub(
	"Virginia opossum",
	"Opossum",
	su13$CommonName
)
su13$CommonName <- gsub(
	"Striped Skunk",
	"Skunk",
	su13$CommonName
)
su13$CommonName <- gsub(
	"Red fox",
	"Redfox",
	su13$CommonName
)

su13 <- su13[-grep("cottontail", su13$CommonName),]

# bring in fall 2013

oc13 <- read.csv(
	"./data/fall_13.csv",
	stringsAsFactors = FALSE
)

oc13$Season <- "FA13"

oc13 <- oc13[,c(1,2,5,3,4)]
colnames(oc13) <- colnames(test3)

# combine all of them

alldat <- bind_rows(test3, su13, oc13)

# Reduce down to only sites in Chicago

palocs <- read.csv(
	"./data/PA_sites.csv",
	stringsAsFactors = FALSE
)

palocs <- sf::st_as_sf(
	palocs,
	coords = c("Easting", "Northing"),
	crs = st_crs(32616)
)

st_crs(palocs) <- st_crs(32616)



library(sf)

county <- sf::st_read(
	"../../GIS/chicago_boundary",
	layer = "geo_export_fbe170bd-acbb-4957-9bc5-7377a2993491"
)


county <- sf::st_transform(
	county,
	crs = st_crs(palocs)
)

ack <- st_buffer(county["name"], 0)

plot(county["name"])
plot(st_geometry(palocs), pch = 15)
plot(ack["name"], add = TRUE)
points(palocs["StationID"])

yo <- st_intersection(palocs, ack)
