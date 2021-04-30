library(uwinutils)
library(lubridate)
connect2db()

# bring in the data

tmp_qry <- "
SELECT sp.CommonName, cl.locationAbbr, ph.photoDateTime, ph.photoName, sa.defaultTimeZone FROM DetectionSpecies ds
INNER JOIN Detections de ON de.detectionID = ds.detectionID
INNER JOIN Photos ph ON ph.photoName = de.photoName
INNER JOIN Visits vi ON vi.visitID = ph.visitID
INNER JOIN CameraLocations cl ON cl.locationID = vi.locationID
INNER JOIN StudyAreas sa ON sa.areaID = cl.areaID
INNER JOIN Species sp ON sp.speciesID = ds.speciesID
WHERE sp.commonName IN ('Coyote', 'Virginia opossum', 'Striped Skunk', 'Raccoon', 'Red fox', 'Eastern cottontail rabbit')
AND sa.areaAbbr = 'CHIL'
AND de.valStatID != 3;"

yo <- SELECT(tmp_qry)

yo$photoDateTime <- with_tz(yo$photoDateTime, 'America/Chicago')
yo$date <- as.Date(yo$photoDateTime, format = "%Y-%M-%d")

# Merge sites that need to be merged
to_combine <- read.table("../Transect/sites_to_merge_sp_13.txt", header = TRUE, sep = "\t")

o1 <- to_combine[to_combine$group==1,]
o1$site_no_number <- paste0(o1$site_no_number, "0")

for(i in 1:nrow(o1)){
	yo$locationAbbr[yo$locationAbbr == o1$Site[i]] <- o1$site_no_number[i]
}

# get the detections
yo2 <- yo %>% group_by(CommonName, locationAbbr, date) %>% 
	summarise(count = length(date)) %>% 
	filter(date > ymd('2013-6-30')) %>% 
	filter(date < ymd('2013-8-1'))

# now we need to get the active range of each camera.

tmp_qry <- "
SELECT cl.locationAbbr, ph.photoDateTime, ph.photoName, sa.defaultTimeZone, vi.visitID FROM Photos ph
INNER JOIN Visits vi ON vi.visitID = ph.visitID
INNER JOIN CameraLocations cl ON cl.locationID = vi.locationID
INNER JOIN StudyAreas sa ON sa.areaID = cl.areaID
WHERE sa.areaAbbr = 'CHIL';"

photos <- SELECT(tmp_qry)
photos$photoDateTime <- with_tz(photos$photoDateTime, "America/Chicago")
photos$date <- as.Date(photos$photoDateTime, format = "%Y-%M-%d")

# Merge sites that need to be merged
to_combine <- read.table("../Transect/sites_to_merge_sp_13.txt", header = TRUE, sep = "\t")

o1 <- to_combine[to_combine$group==1,]
o1$site_no_number <- paste0(o1$site_no_number, "0")

for(i in 1:nrow(o1)){
	photos$locationAbbr[photos$locationAbbr == o1$Site[i]] <- o1$site_no_number[i]
}



start_end <- photos %>% group_by(locationAbbr, visitID) %>% 
	filter(date > ymd('2013-6-30')) %>% 
	filter(date < ymd('2013-8-1')) %>% 
	summarise(min = min(date),
						max = max(date))

# construct the observation matrix
dmin <- min(start_end$min)
dmax <- max(start_end$max)

ndays <- length(seq(dmin, dmax, by = "1 day"))
nsite <- dplyr::n_distinct(photos$locationAbbr)

sites <- unique(photos$locationAbbr)
obs_mat <- matrix(0, ncol = ndays, nrow = nsite)

row.names(obs_mat) <- sites

for(i in 1:nrow(start_end)){
	the_row <- which(sites == start_end$locationAbbr[i])
	to_1 <- day(start_end$min[i]):day(start_end$max[i])
	obs_mat[the_row, to_1] <- 1
	rm(to_1)
}

days_sampled <- rowSums(obs_mat)

days_sampled <- data.frame(locationAbbr = names(days_sampled),
													 J = as.numeric(days_sampled))

all_sites <- expand.grid(CommonName = unique(yo2$CommonName),
									Sites = sites,
									stringsAsFactors = FALSE)

all_sites <- all_sites[order(all_sites$CommonName, all_sites$Sites),]


# now get the number of days detected for each season
y <- yo2 %>% group_by(CommonName, locationAbbr) %>% 
	summarise(count = n_distinct(date)) %>% 
	right_join(., all_sites, by = c("CommonName" = "CommonName", "locationAbbr" = "Sites")) %>% 
	inner_join(., days_sampled, by = "locationAbbr")



to_0 <- which(is.na(y$count) & y$J > 0)
y$count[to_0] <- 0

write.csv(y, "./data/summer_2013.csv", row.names = FALSE, quote = FALSE)
# add zeros where we need them



det_mat <- array(NA)

# Read in the summer data
summer <- data.table::fread("./data/summer_2013.csv", skip = 3, data.table = FALSE)

# Merge sites that need to be merged
to_combine <- read.table("../Transect/sites_to_merge_sp_13.txt", header = TRUE, sep = "\t")

o1 <- to_combine[to_combine$group==1,]
o1$site_no_number <- paste0(o1$site_no_number, "0")

for(i in 1:nrow(o1)){
	summer$Site[summer$Site == o1$Site[i]] <- o1$site_no_number[i]
}

library(dplyr)

# collapse detection history
test <- summer %>% dplyr::group_by(Site, Species) %>% 
	dplyr::select(dplyr::starts_with("Day")) %>% 
	summarise_all(sum)


dh_collapse <- function(x){
	na_count <- sum(is.na(x))
	
	
}




pic_rec$SurveyID <- factor(pic_rec$SurveyID)

event_tbl$site <- substr(event_tbl$SurveyID, 1, 8)
event_tbl$SurveyID <- as.character(event_tbl$SurveyID)
for(i in 1:nrow(o1)){
	event_tbl$SurveyID[event_tbl$site == o1$Site[i]] <- 
		paste0(o1$site_no_number[i], "-", 
					 substr(event_tbl$SurveyID[event_tbl$site == o1$Site[i]], 10,13))
}
event_tbl$SurveyID <- factor(event_tbl$SurveyID)

sur_rec$StationID <- as.character(sur_rec$StationID)
for(i in 1:nrow(o1)){
	sur_rec$StationID[sur_rec$StationID == o1$Site[i]] <- o1$site_no_number[i]
}
stations$StationID <- as.character(stations$StationID)
for(i in 1:nrow(o1)){
	stations$StationID[stations$StationID == o1$Site[i]] <- o1$site_no_number[i]
}
stations <- stations[,c("StationID", "Easting", "Northing")]