##############################################
#
#
# A script for the first figure of this paper
#
# M. Fidino
#
##############################################

# load the plotting and utility functions
source("sourcer.R")

packs <- c(
	"raster",
	"mvtnorm",
	"runjags", 
	"rjags",
	"vioplot",
	"scales",
	"coda",
	"sf"
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

# 
x_cov <- 0.3 * gen_mvn(plane, c(0.2,0.8), c(2,2), 0.1) + 
	0.7 * gen_mvn(plane, c(0.2,0.1), c(1,1), 0.4)

# add this covariate to plane raster object.
values(plane) <- plogis(as.numeric(scale(x_cov)) * 2.5)
names(plane) <- 'x'

plot(
	plane$x,
	col = sf.colors(100, alpha = 0.8)
)


# Make a conflict potential figure

w_cov <-  gen_mvn(plane,c(0.1,0.1), c(0.25,1), 0.1) #+
	#0.5 * gen_mvn(plane, c(0.1,0.3), c(0.2,0.3), 0.4)

test <- rep(seq(-4,4, length.out = 1000), each = 1000)
tmp <- blank

values(tmp) <- plogis(as.numeric(scale(w_cov)) * -2.5)
names(tmp) <- 'w'
plot(tmp$w)

plane <- addLayer(plane, tmp)
rm(tmp)
plot(
	plane$w,
	col = sf.colors(100, alpha = 0.8),
	breaks = seq(0,1,0.01)
)

tmp <- blank
values(tmp) <- values(plane$w) * values(plane$x)
names(tmp) <- "xw"

plot(
	tmp,
	col = sf.colors(100, alpha = 0.8),
	breaks = seq(0,1,0.01)
)

# just going to piece this together in inkscape

windows(6,4)
svg("./figures/rough_1.svg", height = 4, width = 6)
par(mfrow = c(1,4), mar = c(1.5,1.5,1.5,1.5))
plot(plane$x,
		 col = sf.colors(100, alpha = 0.8),
		 breaks = seq(0,1,0.01),
		 axes = FALSE,
		 box = FALSE,
		 legend = FALSE,
		 main = ""
)
plot(plane$w,
		 col = sf.colors(100, alpha = 0.8),
		 breaks = seq(0,1,0.01),
		 axes = FALSE,
		 box = FALSE,
		 legend = FALSE,
		 main = ""
)

plot(tmp, col=sf.colors(100, alpha = 0.8),
		 legend.width=1, legend.shrink=0.75,
		 breaks = seq(0,1,0.01),
		 axes = FALSE,
		 box = FALSE,
		 main = "",
axis.args=list(at=seq(0, 1, 0.25),
							 labels=seq(0, 1, 0.25), 
							 cex.axis=1),
legend.args=list(text='Probability', side=4, font=2, line=2.5, cex=1))
dev.off()
