## List of the R-files supplied with the paper: Koshkina et al. "Integrated Species Distribution Models: Combining presence-only data and presence-absence data with imperfect detection"

PBSO-function.r - Fitting PB, SO and Integrated models the data. The data required for the file is stored in the file data.rda and include

	 - s.occupancy - raster with background covariates that effect occupancy
	 - s.detection - raster with background covariates that effect detection
	 - pb.occupancy - matrix with covariates that effect occupancy in the locations of detected presences of the opportunistic survey
	 - pb.detection - matrix with covariates that effect detection in the locations of detected presences of the opportunistic survey
	 - y.so - matrix of detection/non detection of the SO surveys
	 - so.occupancy - matrix with covariates that effect occupancy in the locations of SO survey sites
	 - so.detection - matrix with covariates that effect detection in the locations of SO survey sites in 		each survey


sim-data.r  - Simulated data generated. Simulates covariates _x_ and _w_, PB and SO datasets.  All the data is saved in a file data.rda

functions.r  - Utility and likelihood functions required for PBSO-function.r

