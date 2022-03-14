source("./R/sourcer.R")
packs <- c(
	"lubridate", "raster", "sp", "sf", "runjags", "coda", "mgcv"
)
package_load(packs)


# read in the data

chicago_raster <- readRDS("./data/all_raw_layers.rds")

# get bounds of plot

chicago_scaled <- raster::scale(
	chicago_raster
)

occ_covs <- raster::values(chicago_scaled)

occ_covs <- occ_covs[complete.cases(occ_covs),]

chicago_raster$URB1 <- NA
chicago_raster$URB2 <- NA

val_locs <- which(!is.na(values(chicago_raster$canopy)))

var_pca <- prcomp(
	occ_covs[,c("canopy", "grass", "imperv", "houses")],
	center = FALSE
)
var_pca$x <- var_pca$x * -1
var_pca$rotation <-  var_pca$rotation * -1

values(chicago_raster$URB1)[val_locs] <- var_pca$x[,1]

values(chicago_raster$URB2)[val_locs] <- var_pca$x[,2]

plot(chicago_raster$URB1)


windows(12,12)

svg("./figures/tmp_supp_map.svg", height = 12, width = 12)
m <- matrix(
	1:6,
	ncol =2,
	nrow = 3,
	byrow = TRUE
)
layout(m)
par(mar = c(5.5,5,4,2))

plot(
	chicago_raster$canopy,
	col = sf.colors(10),
	bty = 'l',
	main = "Canopy cover"
)
plot(
	chicago_raster$grass,
	col = sf.colors(10),
	bty = 'l',
	main = "Grass cover"
)
plot(
	chicago_raster$houses,
	col = sf.colors(10),
	bty = 'l',
	main = "Housing density"
)
plot(
	chicago_raster$imperv,
	col = sf.colors(10),
	bty = 'l',
	main = "Impervious cover"
)

plot(
	chicago_raster$URB1,
	col = sf.colors(10),
	bty = 'l',
	main = "URB1"
)

plot(
	chicago_raster$URB2,
	col = sf.colors(10),
	bty = 'l',
	main = "URB2"
)
dev.off()
