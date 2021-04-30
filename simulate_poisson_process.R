##############################################
#
# Simulate and analyze a poisson point process
#
# Written by M. Fidino 1/9/2020 mdy
#
##############################################

# A lot of this code is based off of:
#  CITATION

############################################
# Setting up the area we are sampling
############################################

# load the plotting and utility functions
source("sourcer.R")

packs <- c(
	"raster",
	"mvtnorm",
	"runjags", 
	"rjags",
	"vioplot",
	"scales",
	"coda"
)

# this will load these packages (and download if necessary.)
package_load(packs)

# Whether or not you want to see the plots from each data
#  simulating step.
do_plots <- TRUE

# bounds of the plane we are sampling within
plane_xmin <- -1
plane_xmax <-  1
plane_ymin <- -1
plane_ymax <-  1

# number of pixels in space. 
npixels_x <- 1000
npixels_y <- 1000

# create a raster, currently has no data. We also make a blank raster
#  so that we can easily add other layers.
plane <- blank <- raster(
	ncol = npixels_x,
	nrow = npixels_y,
  xmn = plane_xmin,
	xmx = plane_xmax,
	ymn=plane_ymin,
	ymx=plane_ymax
)

# the x y coordinates of each pixel in the plane
plane_coord <- xyFromCell(
	plane,
	1:ncell(plane)
)

# generate a spatial covariate. gen_mvn from sim_utility.R script.
x_cov <- 0.4 * gen_mvn(plane, c(0.1,0.8), c(1,1), 0.1) + 
	       0.6 * gen_mvn(plane, c(0.7,0.2), c(0.5,0.75), 0.4)

# add this covariate to plane raster object.
values(plane) <- as.numeric(scale(x_cov))
names(plane) <- 'x'

# Species presence across the landscape. gen_process from sim_utility.R script.
sp_pres <- gen_process(
	rast = plane,
	beta = matrix(
		c(6,1, 4.5, -1, 5, 2),
		ncol = 2,
		nrow = 3,
		byrow=TRUE
	),
	my_seed = 4,
	dm = cbind(
		1, values(plane$x)
	)
)

# plot out species presence across landscape. 
#  plot_dist from plot_utility.R script
if(do_plots){
  plot_dist(
  	plane,
  	cov_name = "x",
  	pixel_id = sp_pres$pixel_id
  )
}

#############################
# simulate presence only data
#############################

# Calculating this at the same resolution as the latent Poisson point process
#  We will later aggregate this to be at the same scale as the presence
#  absence data.

# make a detection covariate that influences presonce only detection
w_cov <- 0.5 * gen_mvn(plane,c(0.7,0.8), c(0.2,0.3), 0.2) +
	0.5 * gen_mvn(plane, c(0.1,0.3), c(0.2,0.3), 0.4)

# add this to the raster
temp <- blank
values(temp) <- as.numeric(scale(w_cov))
names(temp) <- "det"
plane <- addLayer(plane, temp)

# detection linear predictor for thinned Poisson process
beta_det <- c(1, 3)

po_data <- sample_po(
	rast = plane, 
	dm = cbind(1, values(plane$det)),
	beta_det,
	pres = sp_pres
)

# plot out just the presence only data
if(do_plots){
  plot_dist(
  	plane,
  	cov_name = "x",
  	pixel_id = po_data$pixel_id
  )
}

# aggregate down to the smaller scale for presence / absence sampling
agg_factor <- 30
agg_plane <- raster::aggregate(
	plane,
	fact = agg_factor
)

# function to aggregate pixels from gen_process and sample_po
agg_po_pixel_id <- agg_pres(
	plane,
	pixel_id = po_data$pixel_id, 
	agg_factor = agg_factor
)

if(do_plots){
  plot_dist(
  	agg_plane,
  	cov_name = "x",
  	pixel_id = agg_po_pixel_id
  )
}

##################################
# generate presence absence data
##################################

# We need to discretize the sample area to be more course (assuming that
#  our sampling method captures a slightly wider area than how we simulated
#  the latent Poisson process). 

agg_plane <- aggregate(
	plane,
	fact = agg_factor
)
agg_loc <- xyFromCell(
	agg_plane,
	1:ncell(agg_plane)
)

# We need to know which cells the species is and is not in this aggregated cell.
agg_pixelid <- agg_pres(
	plane,
	pixel_id = sp_pres$pixel_id,
	agg_factor
)

# sample the presence absence data
pa_data <- sample_pa(
	agg_plane, 
	n = 300,
	visits = 4, 
	pixel_id = agg_pixelid,
	det_prob = 0.3
)

if(do_plots){
  plot_dist(
  	agg_plane,
  	"x",
  	agg_pixelid
  )
  points(
  	agg_loc[pa_data$site_pixel,],
  	pch = 16,
  	col = "red"
  )
}

# fit an occupancy model with the grid based approach

# the number of grid points
G <- ncell(agg_plane)

# model using just the presence absence data
my_data <- list(
	# Total number of grid points
	G = G,
	# The occupancy covariates of the latent state Poisson process
	occ_covs = cbind(1, values(agg_plane$x)),
	# The detection covariates of the presence absence data
	pa_det_covs = matrix(1, ncol = 1, nrow = G),
	# The pixels the presence only data occurred
	pa_pixel = pa_data$y_mat$pixel,
	# The number of days a species was detected per site / season
	#   for the presence absence data
	y_pa = as.matrix(pa_data$y_mat[,-1]),
	# The number of presence absence sites
	npa = nrow(pa_data$y_mat),
	# The cell area, replicated for EACH cell.
	cell_area = rep(prod(res(agg_plane)), G),
	# The number of seasons sampled
	nyear = nrow(sp_pres$beta),
	# the number of latent parameters
	nlatent = ncol(sp_pres$beta),
	# The number of observational parameters
	nobs = 1
)

m1 <- run.jags(
	model = "presence_absence_pp_dynamic.R", 
	data = my_data, 
	n.chains = 4, 
	inits = inits_pa, 
	monitor = c("beta_occ", "beta_pa_det"), 
	adapt = 1000, 
	burnin = 2000, 
	sample = 3000,
	thin = 2,
	method = 'parallel',
	summarise = FALSE
)

# Fit the integrated model

# get the number of presence only points per season
npo_count <- sapply(agg_po_pixel_id, length)

my_data <- list(
	# Total number of grid points
	G = G, 
	# The occupancy covariates of the latent state Poisson process
	occ_covs = cbind(1, values(agg_plane$x)),
	# The detection covariates for the presence only data
	po_det_covs = cbind(1, values(agg_plane$det)),
	# The detection covariates for the presence / absence data
	pa_det_covs = matrix(1, ncol = 1, nrow = G),
	# The pixels that the presence absence data occurred
	pa_pixel = pa_data$y_mat$pixel,
	# The pixels that the presence only data occurred							
	po_pixel = unlist(agg_po_pixel_id),
	# What season each opportunistic data point came from
	opp_year = rep(1:3, times = npo_count),
	# The number of presence only data points per season
	npo = as.numeric(npo_count),
	# The total number of all presence only data points
	all_npo = sum(npo_count),
	# The number of days a species was detected per site / season
	#   for the presence absence data
	y_pa = as.matrix(pa_data$y_mat[,-1]),
	# The number of sites presence absence data was sampled
	npa = nrow(pa_data$y_mat),
	# The log cell area, logged as the parameters are on the log scale
	#   in the model.
	cell_area = log(prod(res(agg_plane))),
	# Ones for the 'ones' trick in JAGS (for coding up likelihood)
	ones = rep(1, sum(npo_count)),
	# A big constant value for 'ones' trick.
	CONSTANT = 10000,
	# Number of latent parameters
	nlatent = 2,
	# Number of observational parameters, presence only
	nobs_po = 2,
	# Number of observational parameters, presence / absence
	nobs_pa = 1,
  # Number of seasons sampled
	nyear = 3
	)
mstart1 <- Sys.time()
m2 <- run.jags(model = "integrated_pp_simulate.R", 
							 data = my_data, 
							 n.chains = 6, 
							 inits = inits, 
							 monitor = c("beta_occ", "beta_po_det", "beta_pa_det"), 
							 adapt = 1000, 
							 burnin = 2000, 
							 sample = 3000,
							 thin = 5,
							 method = 'parallel',
							 summarise = FALSE)
mtime1 <- Sys.time() - mstart1

# try the faster one, saves about 6 minutes in this model run.

mstart2 <- Sys.time()
m3 <- run.jags(model = "integrated_pp_simulate_faster.R", 
							 data = my_data, 
							 n.chains = 6, 
							 inits = inits, 
							 monitor = c("beta_occ", "beta_po_det", "beta_pa_det"), 
							 adapt = 1000, 
							 burnin = 2000, 
							 sample = 3000,
							 thin = 5,
							 method = 'parallel',
							 summarise = FALSE)
mtime2 <- Sys.time() - mstart2

##########################
# summarize the two models
##########################

# get mcmc matrix from PA model
#  Removing the first column because it just signifies which chain the mcmc
#  step came from.
pamm <- as.matrix(as.mcmc.list(m1), chains = TRUE)[,-1]

# calculate median and 95% CI
pamm_ci <- apply(pamm, 2, quantile, probs = c(0.025,0.5, 0.975))

# do the same for the integrated model
intmm <- as.matrix(as.mcmc.list(m2), chains = TRUE)[,-1]
intmm_ci <- apply(intmm, 2, quantile, probs = c(0.025,0.5, 0.975))
windows(4,4)
tiff("model_comparison.tiff", height = 4, width = 4, units = "in",
		 res = 300, compression = "lzw")
par(mar = c(2,4,1,1))
plot(1~1, bty = "l", ylim = c(0, 7), xlim = c(0.5,2.5), type = "n",
		 xlab = "", ylab = "", xaxt = "n", yaxt = "n")

axis(1, at = seq(1, 2, 1), labels = F, tck = -0.025)
mtext(c("Intercept", "Slope"), 1, line = 0.6, at = c(1,2), cex = 1.5)

axis(2, at = seq(0,7,1), labels = F, tck = -0.025)
axis(2, at = seq(0,7,1/2), labels = F, tck = -0.025/2)
mtext(text = seq(0,7,1),2, line = 0.75, las = 1, at = seq(0,7,1))
mtext("Estimate", 2, line = 2, at = mean(c(0,7)), cex = 1.5)

lines(x = c(0,1.5), y = c(6,6), lty = 3, lwd = 3)
my_vioplot(pamm[,1], add = TRUE, wex = 0.3, at = 0.75, 
					 col = scales::alpha("#FF7E00", 0.4))
my_vioplot(intmm[,1], add = TRUE, wex = 0.3, at = 1.25, 
					 col = scales::alpha("#7e7F8B", 0.4))

lines(x = c(1.5,3), y = c(1,1), lty = 6, lwd = 3)

my_vioplot(pamm[,2], add = TRUE, wex = 0.3, at = 1.75, 
					 col = scales::alpha("#FF7E00", 0.4))
my_vioplot(intmm[,2], add = TRUE, wex = 0.3, at = 2.25, 
					 col = scales::alpha("#7e7F8B", 0.4))

legend("topright", legend = c("PA", "PA & PO"), 
			 col = c(scales::alpha("#FF7E00",0.4),scales::alpha("#7e7F8B",0.4)),
			 				lty = 1, lwd = 7 , bty = "n", title = "Model")
dev.off()



