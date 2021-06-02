#######################################
#
# Plotting the model outputs
#
# written by M. Fidino 
#
#######################################

source("sourcer.R")
packs <- c(
	"lubridate", "raster", "sp", "sf", "runjags", "coda", "dplyr")
package_load(packs)

# Read in the three mcmc output files
mfiles <- list.files(
	"./mcmc_output/",
	pattern = "_model.RDS",
	full.names = TRUE
)

# read in the run.jags files
my_mcmc <- lapply(
	mfiles,
	readRDS
)

my_sum <- lapply(
	my_mcmc,
	function(x)
		round(summary(x),2)
)
names(my_sum) <- c("coyote", "opossum", "raccoon")

species <- "coyote"
# just getting the data together so we have the 
#  occupancy covariates and the like.
source("format_data_for_analysis.R")

cov1 <- sin((2 * pi/4) * (1:my_data$nyear ))
cov2 <- cos((2 * pi/4) * (1:my_data$nyear ))
my_data$temp_covs <- matrix(c(cov1, cov2), ncol = 2, nrow = length(cov1))


# Plot out the regression results for baseline occupancy
#  and conflict.
# get beta_occ and beta_po_det
parms <- lapply(
	my_sum,
	function(x){
		data.frame(
			x[grep("^beta_occ|beta_po_det", row.names(x)), 1:3]
		)
	}
)

for(i in 1:length(parms)){
	parms[[i]]$species <- names(parms)[i]
	parms[[i]]$parameter <- row.names(parms[[i]])
	row.names(parms[[i]]) <- 1:nrow(parms[[i]])
}
# and combine them all
parms <- dplyr::bind_rows(
	parms
)
# order them by parameter and then species
parms <- parms[
	order(parms$parameter, parms$species),
	]

# split into small data.frames for each parameter
parms <- split(
	parms,
	factor(parms$parameter)
)

pretty_parms <- c(
	"URB1|URB2", "Income", "Vacancy"
)




{windows(8,4)
par(
	mfrow = c(1,2),
	mar = c(4,1,2,4)
)


# step 1. plot out occupancy results
plot(
	1~1, type = "n", bty = "n", ylim = c(0.5, 3.5),
	xlim = c(-1, 1), xaxt = "n", yaxt = "n",
	xlab = "", ylab = "", xaxs = "i", yaxs = "i"
)
u <- par("usr")
axis(1, seq(-1,1, 0.5), labels = FALSE, tck = -0.03 )
axis(1, seq(-1,1, 0.25), labels = FALSE, tck = -0.03/2 )
axis(4, 1:5, labels = FALSE, tck = -0.03/2)
mtext(text = sprintf("%.0f",seq(-1,1, 1)), 1, 
			line = 0.5, at = seq(-1,1, 1),las = 1, cex = 1)
mtext(text = sprintf("%.1f",seq(-0.5,0.5, 1)), 1, 
			line = 0.5, at = seq(-0.5,0.5, 1),las = 1, cex = 1)
mtext("Occupancy (log scale)", side = 1,
			line = 2.25, cex = 1.3)
par(xpd = NA)
text(labels = pretty_parms, x = rep(1.3,3) , y = rev(1:3)+0.18,
			cex = 0.8, las = 1, pos = 1)
text(x = u[1]+ abs(u[1] * 0.1), y = u[4],
		 labels = "a)")
par(xpd = FALSE)


for(box in c(2.5, 4.5)){
	rect(xleft = u[1], ybottom = box - 1, xright = u[2], ytop = box,
			 col = '#f0ece9', border = NA)
}
abline(v = 0, lty = 2)
abline(v = u[2])

my_cols <- c("#24d5f7ff", "#5ee38bff", "#ffb226ff")
for(i in 1:3){
	tmp_y <- seq(4 - i-0.25,4 - i+0.25, length.out = 3)
	ci_x <- parms[[i]][rev(1:3),c(1,3)]
	arrows(
		x0 = ci_x$Lower95,
		y0 = tmp_y,
		x1 = ci_x$Upper95,
		y1 = tmp_y,
		length = 0.03,
		angle = 90, 
		code = 3,
		col = my_cols,
		lwd = 3
	)
	med_x <- parms[[i]]$Median[rev(1:3)]
	points(
		med_x, tmp_y,
		pch = c(21,22,23), 
		bg = c("#24d5f7ff", "#5ee38bff","#ffb226ff" ),
		cex = 1.5
		)
}

# step 2. plot out the conflict results
par(mar = c(4,0.5,2,4.5))
plot(
	1~1, type = "n", bty = "l", ylim = c(0.5, 3.5),
	xlim = c(-1, 5), xaxt = "n", yaxt = "n",
	xlab = "", ylab = "", xaxs = "i", yaxs = "i"
)

axis(1, seq(-1,5, 1), labels = FALSE, tck = -0.03 )
axis(1, seq(-1,5, 0.5), labels = FALSE, tck = -0.03/2 )
axis(2, 1:3, labels = FALSE, tck = -0.03/2)
mtext(text = sprintf("%.0f",seq(-1,5, 1)), 1, 
			line = 0.5, at = seq(-1,5, 1),las = 1, cex = 1)
mtext("Conflict (logit scale)", side = 1,
			line = 2.25, cex = 1.3)


u <- par("usr")
for(box in c(2.5, 4.5)){
	rect(xleft = u[1], ybottom = box - 1, xright = u[2], ytop = box,
			 col = '#f0ece9', border = NA)
}
abline(v = 0, lty = 2)
abline(v = u[1])

my_cols <- c("#24d5f7ff", "#5ee38bff", "#ffb226ff")
for(i in 1:3){
	tmp_y <- seq(4 - i-0.25,4 - i+0.25, length.out = 3)
	ci_x <- parms[[i+3]][rev(1:3),c(1,3)]
	arrows(
		x0 = ci_x$Lower95,
		y0 = tmp_y,
		x1 = ci_x$Upper95,
		y1 = tmp_y,
		length = 0.03,
		angle = 90, 
		code = 3,
		col = my_cols,
		lwd = 3
	)
	med_x <- parms[[i+3]]$Median[rev(1:3)]
	points(
		med_x, tmp_y,
		pch = c(21,22,23), 
		bg = c("#24d5f7ff", "#5ee38bff","#ffb226ff" ),
		cex = 1.5
	)
}
par(xpd = NA)
text(x = u[1]+ abs(u[1] * 0.1), y = u[4],
		 labels = "b)")
legend(x = 4.9, y =2.4,
			 legend = rev(c("coyote", "opossum", "raccoon")),
			 pch = rev(21:23), pt.bg = rev(my_cols), pt.cex = 1.5,
			 horiz = FALSE, cex = 0.9, xjust = 0, bty = "n")
}


# plot out the fourier series stuff


mm <- lapply(
	my_mcmc,
	function(x) as.matrix(as.mcmc.list(x))
)

# get model intercept and fourier stuff

four <- lapply(
	mm,
	function(x) x[,grep("po_mu|temp_occ", colnames(x))]
)

names(four) <- names(my_sum)
# do some predictions
mpred <- vector("list", 3)
for(i in 1:3){
mpred[[i]] <- four[[i]] %*% t(cbind(1, my_data$temp_covs)) #+ my_data$cell_area
mpred[[i]] <- plogis(mpred[[i]])#1 - exp(-exp(mpred[[i]]) )
mpred[[i]] <- apply(mpred[[i]], 2, quantile, c(0.025,0.5,0.975))
}

windows(4,8)
par(mfrow = c(3,1))
par(mar = c(4,6,1,1))

xlabel <- rep(c("WI", "SP", "SU", "FA"), 3)

sublabel = c("a) coyote", "b) V. opossum", "c) raccoon" )

my_cols <- c("#24d5f7ff", "#5ee38bff","#ffb226ff" )
my_range <- c(0.4,1)
for(s in 1:3){
	species <- names(four)[s]
	source("format_data_for_analysis.R")
	tmp <- mpred[[s]]
	plot(
		1~1, type = "n", bty = "n", ylim = my_range,
		xlim = c(1, 12), xaxt = "n", yaxt = "n",
		xlab = "", ylab = ""
	)
	u <- par("usr")
	text(x = u[1]+ (u[1]*0.2), y = u[4] - (u[4]*0.1),
			 labels = sublabel[s], pos = 4, cex = 1.5)
	axis(1, seq(1,12, 1), labels = FALSE, tck = -0.03 )
	#axis(1, seq(-1,1, 0.25), labels = FALSE, tck = -0.03/2 )
	axis(2, seq(0.4,1, 0.1), labels = FALSE, tck = -0.03 )
  axis(2, seq(0.4,1, 0.1/2), labels = FALSE, tck = -0.03/2 )
	axis(2, 1:5, labels = FALSE, tck = -0.03/2)
	mtext(text = xlabel, 1, 
				line = 1, at = seq(1,12, 1),las = 1, cex = 1)
	mtext(text = sprintf("%.1f",seq(0.4,1, .2)), 2, 
				line = 0.8, at = seq(0.4,1, .2),las = 1, cex = 1)
	mtext("Season", side = 1,
				line = 2.8, cex = 1.3)
	mtext("Occupancy", side = 2,
				line = 3.8, cex = 1.3)
	x1 <- 1:12
	x2 <- rev(x1)
	y1 <- tmp[1,]
	y2 <- rev(tmp[3,])
  polygon(
  	x = c(x1, x2),
  	y = c(y1,y2),
  	col = scales::alpha(my_cols[s], 0.5),
  	border = NA
  )
  lines(tmp[2,] ~ c(1:12), lwd = 3)
  npo_data <- (my_data$npo +
  						 	sum(my_data$y_pa[,s]>0, na.rm = TRUE)) / my_data$G
  
  lines(npo_data ~ c(1:12))
  
  points(npo_data ~ c(1:12), pch = 19, col = "white", cex = 2)
  points(npo_data ~ c(1:12))
	}


plot(1~1)
plot(1~1)

mc <- as.matrix(as.mcmc.list(m1), chains = TRUE)
# summarise it
mus <- t(apply(mc[,-1], 2, quantile, probs = c(0.025,0.5,0.975)))

mus_r <- round(mus, 2)

ocovs <- colnames(occ_covs)

windows(6,6)
par(mar = c(4,8,1,1))
plot(x = mus_r[1:5,2], y = rev(1:5), pch = 19, xlim = c(-1,1.5),
		 bty = "l", cex = 2, yaxt = "n", ylab = "",
		 xlab = "Parameter estimate (Occupancy)")
for(i in 1:5){
	y <- 6-i
	lines(x = mus_r[i,-2], y = rep(y,2), lwd = 3)
}
axis(2, at = 1:5, labels = rev(ocovs), las = 2)
abline(v = 0, lty = 2)


windows(6,6)
par(mar = c(4,8,1,1))
plot(x = mus_r[8:12,2], y = rev(1:5), pch = 19, xlim = c(-2,1.5),
		 bty = "l", cex = 2, yaxt = "n", ylab = "",
		 xlab = "Parameter estimate (Conflict)")
for(i in 1:5){
	y <- 6-i
	lines(x = mus_r[i+7,-2], y = rep(y,2), lwd = 3)
}
axis(2, at = 1:5, labels = rev(ocovs), las = 2)
abline(v = 0, lty = 2)
#mus <- mus[,2]
windows(30,30)
par(mfrow = c(4,3))
b0s <- mus[grep("psi_mu", row.names(mus)),2]
b1s <- mus[grep("po_mu", row.names(mus)),2]
b2s <- mus[grep("pa_mu", row.names(mus)),2]

windows(30,70)
par(mfrow = c(4,1))

testing <- exp(b0s + ( my_data$occ_covs %*%  mus[1:5,2] ) +my_data$cell_area)
testing2 <- plogis(
	b1s + (my_data$occ_covs %*% mus[8:12,2])
)

testing3 <- plogis(
	b2s + (my_data$pa_det_covs %*% mus[6:7, 2])
	
)


tprob <- 1 - exp(-testing )

longshot <- chicago_raster
longshot$species <- NA
longshot$species[!is.na(values(longshot$id)) ] <- tprob
longshot$conflict <- NA
longshot$conflict[!is.na(values(longshot$id)) ] <- testing2
longshot$species_conflict <- NA
longshot$species_conflict[!is.na(values(longshot$id))] <- tprob * testing2
longshot$species_noconflict <- NA
longshot$species_noconflict[!is.na(values(longshot$id))] <-  tprob * (1 -testing2)
longshot$detection <- NA
longshot$detection[!is.na(values(longshot$id))] <- testing3

longshot$urb1 <- NA
longshot$urb1[!is.na(values(longshot$id))] <- occ_covs[,1]
plot(longshot[["urb1"]])

longshot$urb2 <- NA
longshot$urb2[!is.na(values(longshot$id))] <- occ_covs[,2]
plot(longshot[["urb2"]])



longshot$oo <- NA
count_p <- table(my_data$po_pixel)
longshot$oo[!is.na(values(longshot$id))][as.numeric(names(count_p))] <- as.numeric(count_p)


hm <- as(longshot, 'SpatialPolygonsDataFrame')
p <- st_as_sf(hm, as_points = FALSE, merge = TRUE)


plot(p["species"],
		 main = "coyote occupancy",
		 border = NA,
		 #breaks = seq(0,1,0.1)
)
plot(
	p["conflict"], main = "conflict",
	border = NA,
	breaks = seq(0,1,0.1),
	reset = FALSE
)
ack <- rasterToPoints(longshot$oo)

plot(
	p["species_conflict"], main = "Conflict hotspots",
	border = NA,
	#breaks = seq(0,1,0.1),
	reset = FALSE
)


windows(6,6)
plot(
	p["species_noconflict"],
	main = "opossum non-conflict hotspots",
	border = NA,
	#breaks = seq(0,1,0.1),
	reset = FALSE
)


points(ack, pch = 19, col = scales::alpha("black", 0.5))
plot(longshot[["rcon"]], main = "raccoon occupancy * conflict")
plot(longshot[["rncon"]], main = "raccoon non-conflict")

plot(p)

#points(ack)

}

plot(longshot[["oo"]], add = TRUE)

my_prior <- rlogis(1e6, 0,1)
my_prior <- log(rgamma(1e6, 1,1))

my_prior <- my_prior + 12.42922
hist(my_prior)
hist(1 / (1 + exp(-my_prior)))
hist(1 - exp(-exp(my_prior)))

sink("tmp_mod.R")
cat("
model{
y ~ dgamma(1,1)
k <- log(y)
j <- 1 - exp(-exp(y))
}		
", fill = TRUE)
sink()

kk <- run.jags("tmp_mod.R", monitor = c("y", "k", "j"))

ka <- as.matrix(as.mcmc.list(kk))