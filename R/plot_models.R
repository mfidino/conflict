#######################################
#
# Plotting the model outputs
#
# written by M. Fidino 
#
#######################################

source("./R/sourcer.R")
packs <- c(
	"lubridate", "raster", "sp", "sf", 
	"runjags", "coda", "dplyr", "mgcv")
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


# just getting the data together so we have the 
#  occupancy covariates and the like.
species <- "raccoon"
source("./R/format_data_for_analysis.R")

# Plot out the regression results for baseline occupancy
#  and conflict (This is for what I think will be
#  figure 3).
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

if(!file.exists("./figures/figure_3.tiff")){
tiff(
	"./figures/figure_3.tiff",
	height = 4,
	width = 8,
	units = "in",
	res = 1200,
	compression = "lzw"
)
par(
	mfrow = c(1,2),
	mar = c(4,1,2,4)
)


# step 1. plot out occupancy results
plot(
	1~1, type = "n", bty = "n", ylim = c(0.5, 3.5),
	xlim = c(-0.5, 0.5), xaxt = "n", yaxt = "n",
	xlab = "", ylab = "", xaxs = "i", yaxs = "i"
)
u <- par("usr")
axis(1, seq(-0.5,0.5, 0.5), labels = FALSE, tck = -0.03 )
axis(1, seq(-0.5,0.5, 0.25), labels = FALSE, tck = -0.03/2 )
axis(4, 1:5, labels = FALSE, tck = -0.03/2)
mtext(text = sprintf("%.0f",seq(0,0, 1)), 1, 
			line = 0.5, at = seq(0,0, 1),las = 1, cex = 1)
mtext(text = sprintf("%.1f",seq(-0.5,0.5, 1)), 1, 
			line = 0.5, at = seq(-0.5,0.5, 1),las = 1, cex = 1)
mtext("Occupancy (cloglog scale)", side = 1,
			line = 2.25, cex = 1.3)
par(xpd = NA)
text(labels = pretty_parms, x = rep(0.65,3) , y = rev(1:3)+0.18,
			cex = 0.8, las = 1, pos = 1)
text(x = u[1]+ abs(u[1] * 0.1), y = u[4],
		 labels = "A)")
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
	xlim = c(-1, 4), xaxt = "n", yaxt = "n",
	xlab = "", ylab = "", xaxs = "i", yaxs = "i"
)

axis(1, seq(-1,4, 1), labels = FALSE, tck = -0.03 )
axis(1, seq(-1,4, 0.5), labels = FALSE, tck = -0.03/2 )
axis(2, 1:3, labels = FALSE, tck = -0.03/2)
mtext(text = sprintf("%.0f",seq(-1,4, 1)), 1, 
			line = 0.5, at = seq(-1,4, 1),las = 1, cex = 1)
mtext("Conflict potential (logit scale)", side = 1,
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
text(x = u[1]+ abs(u[1] * 0.3), y = u[4],
		 labels = "B)")
legend(x = 4.0, y =2.4,
			 legend = c("coyote", "opossum", "raccoon"),
			 pch = rev(21:23), pt.bg = rev(my_cols), pt.cex = 1.5,
			 horiz = FALSE, cex = 0.9, xjust = 0, bty = "n")
dev.off()
}


# Plot out GAM stuff. This will be for supplemental material.
mm <- lapply(
	my_mcmc,
	function(x) as.matrix(as.mcmc.list(x))
)

gamb <- lapply(
	mm,
	function(x) x[,grep("^b\\[", colnames(x))]
)

# make gamb an array so we can get average across time.
# get median of each estimate
gamb_med <- lapply(
	gamb,
	function(x) apply(x, 2, median))
	
# get beta occ stuff
beta_occ <- lapply(
	mm,
	function(x) x[,grep("^beta_occ\\[", colnames(x))]
)
beta_occ_med <- lapply(
	beta_occ,
	function(x) apply(x, 2, median)
)
	

names(gamb) <- names(my_sum)
# do some predictions
mpred <- vector("list", 3)
for(i in 1:3){
	tmp <- t(matrix(gamb_med[[i]], nrow = 10, ncol = 12))
	tmp2 <- matrix(beta_occ_med[[i]], nrow = 3, ncol = 12 )
	mpred[[i]] <- 1 - exp(-exp(
		tmp %*% t(my_data$X) + my_data$cell_area +
		t(tmp2) %*% t(my_data$occ_covs)
		))
}

# plot it out
species_raster <- chicago_raster
species_raster$species <- NA
#species_raster$species[!is.na(values(species_raster$id)) ] <- mpred[[1]][1,]

utm_crs <- 26916
pseason <- expand.grid(
	season  = c("Winter", "Spring", "Summer", "Fall"),
	year = 2011:2013
)
# make it a vector instead
pseason <- apply(pseason, 1, paste, collapse = " ")

# The labels of sub-plots	
labs <- matrix(NA, 3,4)
labs[1:12] <- paste0(LETTERS[1:12],")")
for(sp in 1:3){
	svg(
		paste0("./figures/supl_",names(my_sum)[sp],".svg"),
		height = 9,
		width = 9
	)
	par(mfrow = c(4,3), mar = c(1,2.5,2.5,1), xpd = NA)
  for(i in 1:12){
  	# Baseline occupancy
  	species_raster$species[!is.na(values(species_raster$id)) ] <- mpred[[sp]][i,] 
  	r <- species_raster$species > -Inf
  	# convert to polygons
  	edge_plot <- raster::rasterToPolygons(r, dissolve=TRUE, digits = 4)
  
  	plot(
  		species_raster$species,
  		col = sf.colors(10, alpha = 0.8),
  		axes = FALSE,
  		box = FALSE,
  		breaks = seq(0,1,0.1),
  		legend = FALSE,
  		main = pseason[i]
  	)
  	plot(edge_plot, add = TRUE, lwd = 2)
  	u <- par("usr")
  	par(xpd = NA)
  	text(x = u[1]+15000,y = u[4] + 3500, labels = labs[i] , cex = 1.5)
  	if(i == 1){
  		prettymapr::addnortharrow(pos = "topright", padin = c(0.3,0.1), scale = 0.5)
  		prettymapr::addscalebar(plotepsg = utm_crs, style = "ticks",
  								padin =c(0.35, 0), lwd = 2, label.cex = 1.5)
  	}
  }
	dev.off()
}

# get the average across seasons, for figure 2.
# spatial smoothing terms
gamb <- lapply(
	mm,
	function(x) x[,grep("^b\\[", colnames(x))]
)
# Calculate the average across seasons for each coefficient
gam_av <- gamb
av_oc <- diag(3)
for(i in 1:3){
	tmp <- array(
		gamb[[i]],
		dim = c(nrow(gamb[[i]]), 10, 12)
	)
	gam_av[[i]] <- apply(apply(tmp, c(1,2), mean),2,median)
	av_oc[i,] <- quantile(tmp[,1,] + my_data$cell_area, probs =c(0.025,0.5,0.975))
}

# increases
a1 <- diag(3)
for(i in 1:3){
tmp <- array(
		gamb[[i]],
		dim = c(nrow(gamb[[i]]), 10, 12)
	)
a1[i,] <- quantile(
	tmp[,1,] + beta_occ[[i]][,1] + my_data$cell_area,
	probs = c(0.025,0.5,0.975)
) 
}

round(plogis(a1),2)
# report average occupancy
plogis(av_oc)
# model predictions for average
mpred_mu <- mpred

po <- lapply(
	mm,
	function(x) x[,grep("po_mu|^beta_po_det\\[", colnames(x))]
)

po_med <- lapply(
	po,
	function(x) apply(x, 2, median)
)
cpred <- mpred

# get average across seasons
for(i in 1:3){
	#tmp <- t(matrix(gamb_med[[i]], nrow = 10, ncol = 12))
#	tmp <- apply(tmp, 2, mean)
	mpred_mu[[i]] <- 1 - exp(-exp(
		gam_av[[i]] %*% t(my_data$X) + my_data$cell_area +
			beta_occ_med[[i]] %*% t(my_data$occ_covs)
	))
	cpred[[i]] <- plogis(
		po_med[[i]] %*% t(cbind(occ_covs[,-1], 1))
	)
}
utm_crs <- 26916

{svg("./figures/figure2tmp.svg", height = 9, width = 9)
par(mfrow = c(4,3), mar = c(1,2.5,2.5,1), xpd = NA)

labs <- matrix(NA, 3,4)
labs[1:12] <- paste0(LETTERS[1:12],")")
for(i in 1:3){
	# Baseline occupancy
	species_raster$species[!is.na(values(species_raster$id)) ] <- mpred_mu[[i]] 
	r <- species_raster$species > -Inf
	# convert to polygons
	edge_plot <- raster::rasterToPolygons(r, dissolve=TRUE, digits = 4)

	plot(
		species_raster$species,
		col = sf.colors(10, alpha = 0.8),
		axes = FALSE,
		box = FALSE,
		breaks = seq(0,1,0.1),
		legend = FALSE,
		main = pnames[i]
	)
	plot(edge_plot, add = TRUE, lwd = 2)
  u <- par("usr")

  par(xpd = NA)
	text(x = u[1]+15000,y = u[4] + 3500, labels = labs[i,1] , cex = 1.5)
	if(i == 1){
		text(
			x = u[1]+2000,
			y = mean(u[3:4]),
			labels = expression(paste("Pr(",Psi,")")),
			srt = 90,
			cex = 2
		)
	}
	if(i == 1){
	  addnortharrow(pos = "topright", padin = c(0.5,0.1), scale = 0.5)
		addscalebar(plotepsg = utm_crs, style = "ticks",
								padin =c(0.4, -0.1), lwd = 2, label.cex = 1.5)
	}
}
for(i in 1:3){
	# Baseline occupancy
	species_raster$species[!is.na(values(species_raster$id)) ] <- cpred[[i]] 
	r <- species_raster$species > -Inf
	# convert to polygons
	edge_plot <- raster::rasterToPolygons(r, dissolve=TRUE, digits = 4)
	
	plot(
		species_raster$species,
		col = sf.colors(10, alpha = 0.8),
		axes = FALSE,
		box = FALSE,
		breaks = seq(0,1,0.1),
		legend = FALSE,
		main = ""
	)
	plot(edge_plot, add = TRUE, lwd = 2)
	u <- par("usr")

	par(xpd = NA)
	text(x = u[1]+15000,y = u[4] + 3500, labels = labs[i,2] , cex = 1.5)
	if(i == 1){
		text(
			x = u[1]+2000,
			y = mean(u[3:4]),
			labels = expression(paste("Pr(",eta,")")),
			srt = 90,
			cex = 2
		)
	}
}
for(i in 1:3){
	# Baseline occupancy
	species_raster$species[!is.na(values(species_raster$id)) ] <- cpred[[i]] * mpred_mu[[i]] 
	r <- species_raster$species > -Inf
	# convert to polygons
	edge_plot <- raster::rasterToPolygons(r, dissolve=TRUE, digits = 4)

	
	plot(
		species_raster$species,
		col = sf.colors(10, alpha = 0.8),
		axes = FALSE,
		box = FALSE,
		breaks = seq(0,1,0.1),
		legend = FALSE,
		main = ""
	)
	plot(edge_plot, add = TRUE, lwd = 2)
	u <- par("usr")

	par(xpd = NA)
	text(x = u[1]+15000,y = u[4] + 3500, labels = labs[i,3] , cex = 1.5)
	if(i == 1){
		text(
			x = u[1]+2000,
			y = mean(u[3:4]),
			labels = expression(paste("Pr(",Psi,")xPr(",eta,")")),
			srt = 90,
			cex = 2
		)
	}
}
for(i in 1:3){
	# Baseline occupancy
	species_raster$species[!is.na(values(species_raster$id)) ] <- (1 - cpred[[i]]) * mpred_mu[[i]] 
	r <- species_raster$species > -Inf
	# convert to polygons
	edge_plot <- raster::rasterToPolygons(r, dissolve=TRUE, digits = 4)

	
	plot(
		species_raster$species,
		col = sf.colors(10, alpha = 0.8),
		axes = FALSE,
		box = FALSE,
		breaks = seq(0,1,0.1),
		legend = FALSE,
		main = ""
	)
	plot(edge_plot, add = TRUE, lwd = 2)
	u <- par("usr")
	par(xpd = NA)
	text(x = u[1]+15000,y = u[4] + 3500, labels = labs[i,4] , cex = 1.5)
	if(i == 1){
		text(
			x = u[1]+2000,
			y = mean(u[3:4]),
			labels = expression(paste("Pr(",Psi,")x(1-Pr(",eta,"))")),
			srt = 90,
			cex = 2
		)
	}
	
}
dev.off()
}

# plot out conditional conflict

sp_seq <- rep(1:3, each = 3)
col_seq <- rep(c(2:4), 3)
my_grads <- matrix(
	NA,
	ncol = 3,
	nrow = 500
)
my_grads[,1] <- seq(-2,2, length.out = nrow(my_grads))
my_grads[,2] <- seq(0, 150000, length.out = nrow(my_grads))
my_grads[,3] <- seq(0, 50, length.out = nrow(my_grads))
# get scaled ones for prediction
fp <- my_grads
fp[,2] <- (fp[,2] - mean(values(chicago_raster$income), na.rm = TRUE))/
	sd(values(chicago_raster$income), na.rm = TRUE)
fp[,3] <- (fp[,3] - mean(values(chicago_raster$vacancy), na.rm = TRUE))/
	sd(values(chicago_raster$vacancy), na.rm = TRUE)
my_cols <- rev(c("#24d5f7ff", "#5ee38bff", "#ffb226ff"))


{
	tiff("./figures/figure_4.tiff", height = 4, width = 5,
			 units = "in", res = 1200, compression = "lzw")

m2 <- matrix(
	c(
	  10,rep(1,3),rep(c(4,7),each = 3),
	  10,rep(1,3),rep(c(4,7),each = 3),
	  10,rep(2,3),rep(c(5,8),each = 3),
	  10,rep(2,3),rep(c(5,8),each = 3),
	  10,rep(3,3),rep(c(6,9),each = 3),
	  10,rep(3,3),rep(c(6,9),each = 3),
	  rep(10,10)
	 ),
	  ncol = 10,
	  nrow = 7,
	  byrow = TRUE
)

# flip it a bit to plot downwards
layout(m2)
par(mar = c(1.5,1.5,1.5,1.5))
my_axis <- list(
	seq(-2,2,2),
	seq(0,15,5),
	seq(0,5,5)
)
sp_seq <- rep(1:3, 3)
par_seq <- rep(2:4, each = 3)
col_seq <-rep(2:4,3)
my_labs <- c(paste0(LETTERS[1:9],")"))

for(i in 1:9){
	#if(i %in% c(6,7)){
	#	par(mar = c(2,3,2,1))
	#}
	#if(i %in% c(8:9)){
	#	par(mar = c(2,1,2,3))
	#}
	#if(i %in% c(2,3)){
	#	next
	#}
	
	# get the posterior
	tmp_mat <- mm[[sp_seq[i]]][,
										 c(
										 	"po_mu",
										 	"beta_po_det[1]",
										 	"beta_po_det[2]",
										 	"beta_po_det[3]"
										 	)]
	# get predictions
	to_plot <- tmp_mat[,c(1,par_seq[i])] %*% t(cbind(1,fp[,par_seq[i]-1]))
	to_plot <- apply(to_plot,2, quantile, probs = c(0.025,0.5,0.975))

	plot(1~1, type = "n", xaxt = "n", yaxt = "n", ylim = c(0,1),
			 xlim = range(my_grads[,par_seq[i]-1]),
			 xaxs = "i", yaxs="i", bty = "l")
	u <- par("usr")
	to_mult <- c(1,10000,10)[par_seq[i]-1]
	axis(1, my_axis[[par_seq[i]-1]]*to_mult, labels = FALSE, tck = -0.035)
	tmp <- range(my_axis[[par_seq[i]-1]])
	ou <- c(5,7,5)[par_seq[i]-1]

	axis(1, seq(tmp[1]*to_mult, tmp[2]*to_mult, length.out = ou), labels = FALSE, tck = -0.035/2)
	axis(2, seq(0,1,0.25), tck= -0.035, labels = FALSE)
	axis(2, seq(0,1,0.25/2), tck= -0.035/2, labels = FALSE)
	if(i %in% c(7:9)){
		#stop("add mtext")
	}
	if(i == 2){
		mtext("Conditional conflict potential", 2, line = 3.5)
	}
	if(i == 3){
		mtext("URB2", 1, line = 2.75,  cex = 1)
	}
	if(i == 6){
		mtext("Income", 1, line = 2.75, cex = 1)
		mtext(expression("("*"10K"~km^-2~")"), 1, line = 4.75, cex = 1)
	}
	if(i == 9){
		mtext("Vacancy", 1, line = 2.75, cex = 1)
		mtext(expression("("*calls~km^-2~")"), 1, line = 4.75,cex = 1)
	}
	if(i %in% c(1,2,3)){
		mtext(sprintf("%.1f", c(0,0.5,1)), at = c(0,0.5,1),side = 2,las=1,
					line = 0.75)
	}
	if(i == 3){
		mtext(sprintf("%.0f", c(-2,0,2)), at = c(-2,0,2), side = 1, 
					line = 1
		)
	}
	if(i ==6){
		mtext(
			sprintf("%.0f", my_axis[[par_seq[i]-1]]),
			at = my_axis[[par_seq[i]-1]] *to_mult,
			side = 1, 
					line = 1
		)
	}
	if(i ==9){
		mtext(
			sprintf("%.0f", c(0, 25, 50)),
			at = c(0, 25, 50),
			side = 1, 
			line = 1
		)
	}
	x1 <- my_grads[,par_seq[i]-1]
	x2 <- rev(x1)
	y1 <- to_plot[1,]
	y2 <- rev(to_plot[3,])
	polygon(
		c(x1,x2),plogis(c(y1,y2)),
		border = NA,
		col = scales::alpha(my_cols[sp_seq[i]], 0.5)
	)
	lines(
		x = my_grads[,par_seq[i]-1],
		y = plogis(to_plot[2,]),
		lty = c(1,4,3)[sp_seq[i]],
		col = my_cols[sp_seq[i]],
		lwd = 3
	)

	text(x = u[1], y = u[4] - 0.08, my_labs[i], pos = 4 )
	if(i == 7){
		legend("bottomright", c("coyote", "opossum", "raccoon"),
					 lty = c(1,4,3), col = my_cols, lwd = 4, bty = "n",
					 cex =1, seg.len = 3)
	}
}
}
dev.off()

