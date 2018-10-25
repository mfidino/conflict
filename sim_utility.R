##################################
#
# Utility functions for simulation
#
# Written by M. Fidino
# 
##################################

# gen_mvn: 
#  simulate covariate values over cells in a raster
#
# rast = raster object of plane
#
# mu = numeric vector of length 2. Each mu is proportionally where you want
#  the mean value to be. (0,0) is the bottom left, (1,1) is top right.
# 
# sigma = Variance of covariate on x and y axes
#
# rho = Correlation of covariate on x and y axes
gen_mvn <- function(rast = NULL, mu = NULL,
										sigma = NULL, rho = NULL){
	# error checking
	if(class(rast) != "RasterLayer"){
	stop("rast must be a raster object")
	}
	if(length(mu) != 2 | !is.numeric(mu) | any(mu > 1) | any(mu < 0)){
		stop("mu must be a numeric vector of length 2.")
	}
	if(length(sigma) != 2 | !is.numeric(sigma) | any(sigma < 0)){
		stop("Sigma must be a non-negative numeric vector of length 2.")
	}
	if(length(rho) != 1 | !is.numeric(rho)| rho > 1 |rho < -1){
		stop("rho must be a numeric scalar between -1 and 1.")
	}
	
	# get bounds of raster
	bounds <- extent(rast)
	
	# input a proportion of where you want mu to be on x and y
	mu_loc <- c(bounds@xmin + mu[1] * (bounds@xmax - bounds@xmin),
							bounds@ymin + mu[2] * (bounds@ymax - bounds@ymin))
	
	Sigma <- diag(c(sigma[1] * abs(bounds@xmax - bounds@xmin),
									sigma[2] * abs(bounds@ymax - bounds@ymin)))
	# fill the off diagonal
	Sigma[2:3] <- rep(rho * prod(diag(Sigma)))
	
	to_return <- dmvnorm(xyFromCell(rast, 1:ncell(rast)), 
											 mean=mu_loc, 
											 sigma=Sigma)
}

# gen_process: 
#  Uses a PPP to place a species throughout a given landscape.
#  returns the given parameter values, seed used to generate distribution,
#  and the cell a species is located.
#
# rast = raster object of plane
#
# beta = conformable vector or matrix to values in rast object.
# 
# my_seed = seed to randomly generate values. For reproducibility.
#
# return_occ = return derived occupancy probability in each cell. 
#  defaults to false.

gen_process <- function(rast = NULL, beta = NULL, my_seed = NULL,
												return_occ = FALSE){
	if(class(rast) != "RasterLayer"){
		stop("rast must be a raster object")
	}
	if(is.null(my_seed)){
		my_seed <- floor(runif(1, -1e4, 1e4))
	}
	set.seed(my_seed)
	
	# log-linear intensity
	lin_pred <- exp(as.numeric(cbind(1, values(rast)) %*% beta ))
	
	# get area of raster
	rast_extent <- extent(rast)
	rast_area <- (rast_extent[2] - rast_extent[1]) * 
		(rast_extent[4] - rast_extent[3])
	
	# generate number of individuals potentially on plane
	my_n <- rpois(1, max(lin_pred) * rast_area)
	
	if(my_n > ncell(rast)){
		stop("\nThe number of cells occupied is greater than the total number of cells. \nEither add more cells or decrease beta values.")
	}
	#  sampling w/o replacement ensures only 1 indiv per pixel
	plane_loc <- sort(sample(1:ncell(plane), size=my_n, replace=FALSE))
	
	# put species on landscape relative to how well each cell is for the species
	sp_there = runif(my_n, 0,1) <= lin_pred[plane_loc]/max(lin_pred)
	
	
	# the pixels that have the species
	pixel_id <- plane_loc[sp_there]
	
	# coordinates that have speices
	sp_coord = plane_coord[pixel_id,]
	
	# calcualte it as an occupancy probability via comp log log link
	#  multiply by area of each cell to approximate lambda
	my_occ <- 1 - exp(-lin_pred * prod(res(plane)))
	
	to_return <- list(beta = beta,
										seed = my_seed,
										pixel_id = pixel_id)
	# return derived occupancy probability for each cell
	if(return_occ){
		to_return <- c(to_return, occ_prob)
	}
	return(to_return)
}


# sample pa:
#  

sample_pa <- function(rast = NULL, n = NULL, 
											det_prob = NULL, visits = NULL,
											sp_pixel = NULL, my_seed = NULL){
	if(is.null(my_seed)){
		my_seed <- floor(runif(1, -1e4, 1e4))
	}
	
	# dimensions of raster
	w <- dim(rast)[1]
	h <- dim(rast)[2]
	
	# very rough and somewaht even spacing of cameras
	my_cams <- floor(seq(1, ncell(plane), length.out = n))
	
	# move half the cameras a little bit so they don't end up in a line
	jiggle_cam <- sort(sample(1:n, floor(n/2), replace = FALSE))
	
	my_cams[jiggle_cam] <- floor( my_cams[jiggle_cam] + median(diff(my_cams)/2))
	
	
	# keep all cams within realistic pixel range
	if(any(my_cams > ncell(rast))){
		my_cams[my_cams > ncell(rast)] <- ncell(rast)
	}
	
	
	sp_pres_det <- sort(sp_pres$pixel_id[which(sp_pres$pixel_id %in% my_cams)])
	
	y_mat <- data.frame(pixel = sort(my_cams), y = 0)
	
	y_mat$y[which(y_mat$pixel %in% sp_pres_det)] <- rbinom(length(sp_pres_det), 
																												 size = visits, 
																												 det_prob)
	
	to_return <- list(y_mat = y_mat,
										visits = visits,
										site_pixel = my_cams,
										det_prob = det_prob,
										seed = my_seed,
										sites_detected = length(sp_pres_det))
	return(to_return)
}