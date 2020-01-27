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

