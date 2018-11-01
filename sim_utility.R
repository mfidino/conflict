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
# dm = design matrix. First column must be 1
#
# beta = conformable vector or matrix to values in dm object.
# 
# my_seed = seed to randomly generate values. For reproducibility.
#
# return_occ = return derived occupancy probability in each cell. 
#  defaults to false.

gen_process <- function(rast = NULL, dm = NULL, beta = NULL, my_seed = NULL,
												return_occ = FALSE){
	if(is.null(my_seed)){
		my_seed <- floor(runif(1, -1e4, 1e4))
	}
	set.seed(my_seed)
	
	# cell area on log scale
	cell_area <- log(prod(res(rast)))
	# coords
	rast_coord <- xyFromCell(rast, 1:ncell(rast))
	
	# log-linear intensity
	lin_pred <- exp(as.numeric(dm %*% beta ) + cell_area)
	
	# get area of raster
	rast_extent <- extent(rast)
	rast_area <- (rast_extent[2] - rast_extent[1]) * 
		(rast_extent[4] - rast_extent[3])
	
	# probability of an individual being in a cell
	my_prob <- 1 - exp(-lin_pred)
	
	# expected population size
	pop_size <- sum(lin_pred)
	
	# add some variability to population size
	pop_size <- rpois(1, pop_size)
	
	# place individuals on the landscape relative to the quality of each
	#  cell.
	pixel_id <- sort(sample(x = 1:ncell(rast), 
													size = pop_size, 
													prob = my_prob, 
													replace = FALSE))
	
	# coordinates that have speices
	sp_coord <- rast_coord[pixel_id,]
	
	to_return <- list(beta = beta,
										seed = my_seed,
										pixel_id = pixel_id,
										lambda = lin_pred)
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
	my_cams <- floor(seq(1, ncell(rast), length.out = n))
	
	# move half the cameras a little bit so they don't end up in a line
	jiggle_cam <- sort(sample(1:n, floor(n/2), replace = FALSE))
	
	my_cams[jiggle_cam] <- floor( my_cams[jiggle_cam] + median(diff(my_cams)/2))
	
	
	# keep all cams within realistic pixel range
	if(any(my_cams > ncell(rast))){
		my_cams[my_cams > ncell(rast)] <- ncell(rast)
	}
	
	
	sp_pres_det <- sort(sp_pixel[which(sp_pixel %in% my_cams)])
	
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


# sample_po: sample presence only data
#  returns a vector of the pixels that a species is detected as 'present' 
#  based off of the 'det' covariate layer in object 'rast' and the 'beta' values.
#  Must also include the true presence of indviduals across the landscape that
#  was calculated via gen_process.
sample_po <- function(rast = NULL, dm = NULL, beta = NULL, 
											pres = NULL, my_seed = NULL ){
	if(is.null(my_seed)){
		my_seed <- floor(runif(1, -1e4, 1e4))
	}
	set.seed(my_seed)
	# cell area in log scale
	cell_area <- log(prod(res(rast)))
	# coords
	rast_coord <- xyFromCell(rast, 1:ncell(rast))
	
	#
	lp_det <- plogis(dm %*% beta)
	lp_det[-sp_pres$pixel_id] <- 0
	
	pixel_id <- rbinom(n = ncell(rast), 1, prob = lp_det)
	pixel_id <- which(pixel_id > 0)
	
	# drop half of the points
	pixel_id <- sample(pixel_id, size = floor(length(pixel_id)/1.2), replace = FALSE)
	
	to_return <- list(pixel_id = pixel_id, seed = my_seed, beta = beta)
	return(to_return)
}


sample_po2 <- function(rast = NULL, dm = NULL, beta = NULL, 
											pres = NULL, my_seed = NULL ){
	if(is.null(my_seed)){
		my_seed <- floor(runif(1, -1e4, 1e4))
	}
	set.seed(my_seed)
	# cell area in log scale
	cell_area <- log(prod(res(rast)))
	# coords
	rast_coord <- xyFromCell(rast, 1:ncell(rast))
	
	#
	sp_pres_lam <- log(sp_pres$lambda)
	lp_det <- dm %*% beta
	lin_pred <- sp_pres_lam + lp_det
	my_prob <- as.numeric(1 - exp(-exp(lin_pred)))
	my_n <- sum(sp_pres$lambda)
	my_n <- rpois(1, my_n)
	
	pixel_id <- sample(1:ncell(rast), my_n, prob = my_prob, replace = FALSE)
	pixel_id <- sort(sample(pixel_id, floor(length(pixel_id) * 0.7), replace = FALSE))
	
	
	to_return <- list(pixel_id = pixel_id, seed = my_seed, beta = beta)
	return(to_return)
}

# initial values for presence absence

inits_pa <- function(chain){
	gen_list <- function(chain = chain){
		list( 
			z = rep(1, G),
			beta_occ = rnorm(2),
			beta_pa_det = rnorm(1),
			.RNG.name = switch(chain,
												 "1" = "base::Wichmann-Hill",
												 "2" = "base::Marsaglia-Multicarry",
												 "3" = "base::Super-Duper",
												 "4" = "base::Mersenne-Twister",
												 "5" = "base::Wichmann-Hill",
												 "6" = "base::Marsaglia-Multicarry",
												 "7" = "base::Super-Duper",
												 "8" = "base::Mersenne-Twister"),
			.RNG.seed = sample(1:1e+06, 1)
		)
	}
	return(switch(chain,           
								"1" = gen_list(chain),
								"2" = gen_list(chain),
								"3" = gen_list(chain),
								"4" = gen_list(chain),
								"5" = gen_list(chain),
								"6" = gen_list(chain),
								"7" = gen_list(chain),
								"8" = gen_list(chain)
	)
	)
}

# initial values for integrated sdm

inits <- function(chain){
	gen_list <- function(chain = chain){
		list( 
			z = rep(1, G),
			beta_occ = rnorm(2),
			beta_pa_det = rnorm(1),
			beta_po_det = rnorm(2),
			.RNG.name = switch(chain,
												 "1" = "base::Wichmann-Hill",
												 "2" = "base::Marsaglia-Multicarry",
												 "3" = "base::Super-Duper",
												 "4" = "base::Mersenne-Twister",
												 "5" = "base::Wichmann-Hill",
												 "6" = "base::Marsaglia-Multicarry",
												 "7" = "base::Super-Duper",
												 "8" = "base::Mersenne-Twister"),
			.RNG.seed = sample(1:1e+06, 1)
		)
	}
	return(switch(chain,           
								"1" = gen_list(chain),
								"2" = gen_list(chain),
								"3" = gen_list(chain),
								"4" = gen_list(chain),
								"5" = gen_list(chain),
								"6" = gen_list(chain),
								"7" = gen_list(chain),
								"8" = gen_list(chain)
	)
	)
}


my_vioplot <- function (x, ..., range = 1.5, h = NULL, ylim = NULL, names = NULL, 
												horizontal = FALSE, col = "magenta", border = "black", lty = 1, 
												lwd = 1, rectCol = "black", colMed = "white", pchMed = 19, 
												at, add = FALSE, wex = 1, drawRect = TRUE) 
{
	datas <- list(x, ...)
	n <- length(datas)
	if (missing(at)) 
		at <- 1:n
	upper <- vector(mode = "numeric", length = n)
	lower <- vector(mode = "numeric", length = n)
	q1 <- vector(mode = "numeric", length = n)
	q3 <- vector(mode = "numeric", length = n)
	med <- vector(mode = "numeric", length = n)
	base <- vector(mode = "list", length = n)
	height <- vector(mode = "list", length = n)
	baserange <- c(Inf, -Inf)
	args <- list(display = "none")
	if (!(is.null(h))) 
		args <- c(args, h = h)
	for (i in 1:n) {
		data <- datas[[i]]
		data.min <- min(data)
		data.max <- max(data)
		q1[i] <- quantile(data, 0.25)
		q3[i] <- quantile(data, 0.75)
		med[i] <- median(data)
		iqd <- q3[i] - q1[i]
		upper[i] <- min(q3[i] + range * iqd, data.max)
		lower[i] <- max(q1[i] - range * iqd, data.min)
		est.xlim <- c(min(lower[i], data.min), max(upper[i], 
																							 data.max))
		smout <- do.call("sm.density", c(list(data, xlim = est.xlim), 
																		 args))
		hscale <- 0.4/max(smout$estimate) * wex
		base[[i]] <- smout$eval.points
		height[[i]] <- smout$estimate * hscale
		t <- range(base[[i]])
		baserange[1] <- min(baserange[1], t[1])
		baserange[2] <- max(baserange[2], t[2])
	}
	if (!add) {
		xlim <- if (n == 1) 
			at + c(-0.5, 0.5)
		else range(at) + min(diff(at))/2 * c(-1, 1)
		if (is.null(ylim)) {
			ylim <- baserange
		}
	}
	if (is.null(names)) {
		label <- 1:n
	}
	else {
		label <- names
	}
	boxwidth <- 0.05 * wex
	if (!add) 
		plot.new()
	if (!horizontal) {
		if (!add) {
			plot.window(xlim = xlim, ylim = ylim)
			axis(2)
			axis(1, at = at, label = label)
		}
		# box()
		for (i in 1:n) {
			polygon(c(at[i] - height[[i]], rev(at[i] + height[[i]])), 
							c(base[[i]], rev(base[[i]])), col = col, border = border, 
							lty = lty, lwd = lwd)
			if (drawRect) {
				lines(at[c(i, i)], c(lower[i], upper[i]), lwd = lwd, 
							lty = lty)
				rect(at[i] - boxwidth/2, q1[i], at[i] + boxwidth/2, 
						 q3[i], col = rectCol)
				points(at[i], med[i], pch = pchMed, col = colMed)
			}
		}
	}
	else {
		if (!add) {
			plot.window(xlim = ylim, ylim = xlim)
			axis(1)
			axis(2, at = at, label = label)
		}
		# box()
		for (i in 1:n) {
			polygon(c(base[[i]], rev(base[[i]])), c(at[i] - height[[i]], 
																							rev(at[i] + height[[i]])), col = col, border = border, 
							lty = lty, lwd = lwd)
			if (drawRect) {
				lines(c(lower[i], upper[i]), at[c(i, i)], lwd = lwd, 
							lty = lty)
				rect(q1[i], at[i] - boxwidth/2, q3[i], at[i] + 
						 	boxwidth/2, col = rectCol)
				points(med[i], at[i], pch = pchMed, col = colMed)
			}
		}
	}
	invisible(list(upper = upper, lower = lower, median = med, 
								 q1 = q1, q3 = q3))
}

