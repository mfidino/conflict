

# for raccoon it looks like chain = 2 should be dropped

#m1 <- as.mcmc.list(m1)[-2]

#saveRDS(m1, "./mcmc_output/raccoon_ranef_burnin_trial_kmres.rds")



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