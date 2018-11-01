# Simulate and analyze a poisson point process



############################################
# Setting up the area we are sampling
############################################

# load the utility functions
source("sim_utility.R")

packs <- c("raster", "mvtnorm", "runjags", "rjags", "vioplot", "scales")

# this will load these packages (and download if necessary.)
package_load(packs)

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
plane <- blank <- raster(ncol = npixels_x, nrow = npixels_y,
								         xmn = plane_xmin, xmx = plane_xmax, 
								         ymn=plane_ymin, ymx=plane_ymax)

# the x y coordinates of each pixel in the plane
plane_coord <- xyFromCell(plane, 1:ncell(plane))

# generate a spatial covariate. gen_mvn from sim_utility.R script.
x_cov <- 0.4 * gen_mvn(plane, c(0.1,0.8), c(1,1), 0.1) + 
	       0.6 * gen_mvn(plane, c(0.7,0.2), c(0.5,0.75), 0.4)

# add this covariate to plane raster object.
values(plane) <- as.numeric(scale(x_cov))
names(plane) <- 'x'

# Species presence across the landscape. gen_process from sim_utility.R script.
sp_pres <- gen_process(rast = plane,beta = c(6,1),my_seed = 4,
												dm = cbind(1, values(plane$x)))

# plot out species presence across landscape
plot(plane$x)
points(plane_coord[sp_pres$pixel_id, ], pch = 16)

#############################
# simulate presence only data
#############################

# Calculating this at the same resolution as the latent Poisson process

# make a detection covariate that influences po detection
w_cov <- 0.5 * gen_mvn(plane, c(0.7,0.8), c(0.2,0.3), 0.2) +
	0.5 * gen_mvn(plane, c(0.1,0.3), c(0.2,0.3), 0.4)

# add this 
temp <- blank
values(temp) <- as.numeric(scale(w_cov))
names(temp) <- "det"
plane <- addLayer(plane, temp)

# detection lin_pred for thinned Poisson process
beta_det <- c(1, 3)

po_data <- sample_po(rast = plane, 
										 dm = cbind(1, values(plane$det)),
										 beta_det, pres = sp_pres)

plot(plane$det)
points(plane_coord[po_data$pixel_id,], pch = 16, col = "red")

# aggregate down to the smaller scale for presence / absence sampling
#  For this simulation we are reducing down to 625 cells.
agg_factor <- 10
agg_plane <- aggregate(plane, fact = agg_factor)

# We need to know which cells the species is and is not in this aggregated cell.
temp <- blank
tmp_vals <- rep(0, ncell(temp))
tmp_vals[po_data$pixel_id] <- 1
values(temp) <- tmp_vals
names(temp) <- "z"

agg_po <- aggregate(temp, fact = agg_factor, fun = sum)

agg_loc <- xyFromCell(agg_plane, 1:ncell(agg_plane))

agg_po_pixel_id <- which(values(agg_po)>0)


plot(agg_plane$x )
points(agg_loc[agg_po_pixel_id,])

##################################
# generate presence absence data
##################################

# We need to discretize the sample area to be more course (assuming that
#  our sampling method captures a slightly wider area than how we simulated
#  the latent Poisson process). 

agg_plane <- aggregate(plane, fact = agg_factor)

# We need to know which cells the species is and is not in this aggregated cell.

temp <- blank
tmp_vals <- rep(0, ncell(temp))
tmp_vals[sp_pres$pixel_id] <- 1
values(temp) <- tmp_vals
names(temp) <- "z"

agg_pres <- aggregate(temp, fact = agg_factor, fun = sum)

agg_loc <- xyFromCell(agg_plane, 1:ncell(agg_plane))

agg_pixelid <- which(values(agg_pres)>0)

plot(agg_plane$x)
points(agg_loc[agg_pixelid,])

pa_data <- sample_pa(agg_plane, n = 300, visits = 4, 
									sp_pixel = agg_pixelid,det_prob = 0.3)
points(agg_loc[pa_data$site_pixel,], pch = 16)

# fit an occupancy model with the grid based approach

G <- ncell(agg_plane)

# model using just the presence absence data

my_data <- list(G = G, occ_covs = cbind(1, values(agg_plane$x)),
								pa_det_covs = matrix(1, ncol = 1, nrow = G),
								pa_pixel = pa_data$y_mat$pixel,
								y_pa = pa_data$y_mat$y,
								cell_area = rep(res(agg_plane)[1]*res(agg_plane)[2], 
																length(values(agg_plane$x))),
								npa = nrow(pa_data$y_mat))



m1 <- run.jags(model = "PP_presence_absence_data.R", 
							 data = my_data, 
							 n.chains = 4, 
							 inits = inits_pa, 
							 monitor = c("beta_occ", "beta_po_det", "beta_pa_det",
							 						"test1", "test2", "test3", "test4"), 
							 adapt = 1000, 
							 burnin = 10000, 
							 sample = 10000,
							 method = 'parallel',
							 summarise = FALSE)

# try the integrated model



my_data <- list(G = G, occ_covs = cbind(1, values(agg_plane$x)),
								po_det_covs = cbind(1, values(agg_plane$det)),
								pa_det_covs = matrix(1, ncol = 1, nrow = G),
								pa_pixel = pa_data$y_mat$pixel,
								po_pixel = agg_po_pixel_id,
								y_pa = pa_data$y_mat$y,
								ones = rep(1, length(agg_po_pixel_id)),
								cell_area = log(prod(res(agg_plane))),
								npa = nrow(pa_data$y_mat),
								npo = length(y_po),
								CONSTANT = 10000)

m2 <- run.jags(model = "integrated_ones.R", 
							 data = my_data, 
							 n.chains = 4, 
							 inits = inits, 
							 monitor = c("beta_occ", "beta_po_det", "beta_pa_det"), 
							 adapt = 1000, 
							 burnin = 10000, 
							 sample = 10000,
							 method = 'parallel',
							 summarise = FALSE)

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



