library(uwinutils)
library(lubridate)
connect2db()

# bring in the OC13 data
oc13_fp <- list.files("./data/OC13/", full.names = TRUE)

oc13 <- lapply(oc13_fp, read.csv, stringsAsFactors = FALSE)

# condense to species data

oc13sp <- oc13[c(1,4,5,7:8)]
names(oc13sp) <- c("Coyote", "Raccoon", "Redfox", "Skunk", "Opossum")

test <- lapply(oc13sp, function(x) rowSums(x[,-1], na.rm = TRUE))


ocdf <- data.frame(
	commonName = rep(names(oc13sp), each = length(oc13sp$Coyote$X)),
	locationAbbr = rep(oc13sp$Coyote$X, times = length(oc13sp)),
	count = unlist(test),
	stringsAsFactors = FALSE
)
	

# calculate J

ocdf$J <- rowSums(oc13[[3]][,-1])

# Merge sites that need to be merged
to_combine <- read.table("../Transect/sites_to_merge_sp_13.txt", header = TRUE, sep = "\t")

o1 <- to_combine[to_combine$group==1,]
o1$site_no_number <- paste0(o1$site_no_number, "0")

for(i in 1:nrow(o1)){
	ocdf$locationAbbr[ocdf$locationAbbr == o1$Site[i]] <- o1$site_no_number[i]
}
ocdf$count[ocdf$J == 0] <- NA

write.csv(ocdf, "./data/fall_13.csv", row.names = FALSE, quote = FALSE)
