###################################
#
# Cleaning up the conflict data
#
# Written by M. Fidino 1/6/2020 mdy
#
###################################

# Note: This script assumes you are in the working directory and that the
#  data sub-folder contains:
#  data/conflict_raw
#  data/conflict_clean

# The goal is to read in the data from the raw conflict folder, clean it, and
#  save the cleaned data to conflict_clean
# load packages
library(dplyr)

# load utility script for helper functions
source("utility_script.R")

# Read in the hw conflicts, update the column headers, and clean the data
raccoon <- read.csv(
	file = "./data/conflicts_raw/raccoon.csv",
	stringsAsFactors = FALSE
) %>% 
	update_foia_columns()

# clean the data


opossum <- read.csv(
	file = "./data/conflicts_raw/opossum.csv",
	stringsAsFactors = FALSE
) %>% 
	update_foia_columns()

# drop the first 3 rows
opossum <- opossum[-c(1:3),]

opossum <- read.csv("./data/conflict_clean/opossum.csv")

drop_uncertain_species(opossum, "opossum")

coyote <- read.csv(
	file = "./data/conflicts_raw/coyote.csv",
	stringsAsFactors = FALSE,
) %>% 
	update_foia_columns()

# drop the first 3 rows, which are blank
coyote <- coyote[-c(1:3),]

# go through and clean descriptions.
drop_uncertain_species(coyote, "coyote")
