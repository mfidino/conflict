# Simulate a spatial covariate

library(raster)
library(mvtnorm)

source("sim_utility.R")

# bounds of the plane we are sampling within
plane_xmin = -1
plane_xmax =  1
plane_ymin = -1
plane_ymax =  1

# area of space
plane_area =  (plane_xmax-plane_xmin)*(plane_ymax-plane_ymin)

# number of pixels in space
npixels_x = 30
npixels_y = 30

# create a raster, currently has no data
plane = raster(ncol=npixels_x, nrow=npixels_y, xmn=plane_xmin, xmx=plane_xmax, 
					 ymn=plane_ymin, ymx=plane_ymax)

# the x y coordinates of each pixel in the plane
plane_coord = xyFromCell(plane, 1:ncell(plane))

# generate a spatial covariate

x_cov <- 0.4 * gen_mvn(plane, c(0.1,0.8), c(1,1), 0.1) + 
	       0.6 * gen_mvn(plane, c(0.7,0.2), c(0.5,0.75), 0.4)

values(plane) = as.numeric(scale(x_cov))
names(plane) = 'x'
plot(plane)


sp_pres <- gen_process(plane, c(2,2), 4)


# plot out species presence across landscape
plot(plane$x)
points(plane_coord[sp_pres$pixel_id, ], pch = 16)



pa_data <- sample_pa(plane,200, 0.4, visits = 4)
# push every other row over by delta x?

plot(plane$x)
points(plane_coord[sp_pres$pixel_id,], pch = 16)
points(plane_coord[as.numeric(my_cams),], col = "brown")
# do camera trapping

my_sites <- as.numeric(my_cams)

# detection probability
det <- 0.4




points(plane_coord[y_mat$pixel[y_mat$y > 0],], pch = 16, col = "blue")
# get presence only data

# make a detection covariate that influences po detection

w_cov <- 0.2 * gen_mvn(s, c(0.2,0.8), c(0.5,0.75), 0.3) + 
	0.8 * gen_mvn(s, c(0.3,0.9), c(0.5,0.3), 0.4)

temp <- s
values(temp) <- as.numeric(scale(w_cov))
names(temp) <- "det"
s <- addLayer(s, temp)



# detection lin_pred
beta_det <- c(-0.4, 2)

lp_det <- as.numeric(beta_det %*% t(cbind(1, values(s$det))))

# which areas do we have a probability of getting po
det_prob_po <- plogis(lp_det)
det_prob_po[-pixel_id] <- 0

po <- rbinom(length(det_prob_po), 1, det_prob_po) 

po_pixel <- c(1:400)[po == 1]
plot(s$x, bty = "n")
points(s.loc[pixel_id,], pch = 16, col = "black")
points(s.loc[as.numeric(my_cams),], col = "brown")
points(s.loc[po_pixel,], pch = 16, col = "red")
points(s.loc[y_mat$pixel[y_mat$y > 0],], pch = 16, col = "blue")

# fit an occupancy model with the grid based approach

G <- ncell(plane)


my_data <- list(G = G, occ_covs = cbind(1, values(plane$x)),
								det_covs = matrix(1, ncol = 1, nrow = G),
								pixel = y_mat$pixel,
								y = y_mat$y,
								cell_area = rep(res(plane)[1]*res(plane)[2], 
																length(values(plane$x))),
								nsite = nrow(y_mat))



library(runjags)
library(rjags)
load.module('glm')

inits <- function(chain){
	gen_list <- function(chain = chain){
		list( 
			z = rep(1, G),
			beta_occ = rnorm(2),
			beta_det = rnorm(1),
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

m1 <- run.jags(model = "PP_presence_absence_data.R", 
							 data = my_data, 
							 n.chains = 4, 
							 inits = inits, 
							 monitor = c("beta_occ", "beta_det", "psi"), 
							 adapt = 1000, 
							 burnin = 5000, 
							 sample = 10000,
							 method = 'parallel',
							 summarise = FALSE)

m2 <- as.matrix(as.mcmc.list(m1), chains = TRUE)
zz <- m2[,5:ncol(m2)]

zz <- apply(zz, 2, quantile, probs = c(0.025,0.5,0.975))
zor <- order(zz[2,])

plot(zz[2,zor], type = 'l', ylim = c(0,1), bty = "l", lwd = 3)
lines(zz[1,zor], lty =3)
lines(zz[3,zor], lty =3)

points(my_occ[zor])
plot(zz[2,] ~ lin_pred, type = 'l')
zz <- colMeans(zz)

parms <- m2[,2:4]

apply(parms, 2, quantile, probs = c(0.025,0.5,0.975))
temp <- plane
temp <- dropLayer(temp, 1)

values(temp) <- my_occ
names(temp) <- "psi2"
plane <- addLayer(plane, temp)

plot(plane$psi2)
points(sp_coord, pch = 16)
# make a subset presence only

po <- sample(pixel_id, 75, replace = FALSE)
points(s.loc[po,], pch = 16, col = "red")
points(s.loc[y_mat$pixel[y_mat$y > 0],], pch = 16, col = "blue")


mu2.x = s.xmin + 0.15*(s.xmax-s.xmin)
mu2.y = s.ymin + 0.80*(s.ymax-s.ymin)
sigma2.x = 0.50*abs(s.xmax-s.xmin)
sigma2.y = 0.25*abs(s.ymax-s.ymin)
rho2.xy = -0.4
mu2 = c(mu2.x, mu2.y)
Sigma2 = matrix(c(sigma2.x^2, rep(rho2.xy*sigma2.x*sigma2.y, 2), sigma2.y^2), ncol=2)

xcov = 0.4 * dmvnorm(s.loc, mean=mu1, sigma=Sigma1) + 0.6 * dmvnorm(s.loc, mean=mu2, sigma=Sigma2)
xcov = (xcov - mean(xcov))/sd(xcov)