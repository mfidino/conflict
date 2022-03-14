source("./R/sourcer.R")
packs <- c(
	"lubridate", "raster", "sp", "sf", 
	"runjags", "coda", "dplyr", "mgcv")
package_load(packs)



# needed to subset the data for validating
#  in the formatting data script.
validation_time <- TRUE

# Read in the three mcmc output files
mfiles <- list.files(
	"./mcmc_output/validation/",
	pattern = "validation_model.RDS",
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

grab <- function(x,y) x[,grep(y, colnames(x))]

# This will just grab the validation data

brier <- matrix(
	NA,
	ncol = 3,
	nrow = nrow(do.call("rbind", my_mcmc[[1]]$mcmc))
)
for(sp in 1:3){
	species <- names(my_sum)[sp]
  source("./R/format_data_for_validation.R")
	
	# generate predictions for the validation data
	mc <- do.call("rbind", my_mcmc[[sp]]$mcmc)
	con_conflict <- vector("list", length = my_data$nyear)
	pb <- txtProgressBar(max = my_data$nyear)
	for(year in 1:my_data$nyear){
		setTxtProgressBar(pb, year)
		# Get latent occupancy across the city
		latent <- my_data$X %*% t(grab(mc, paste0("b\\[\\w\\w?,",year,"\\]"))) +
			my_data$occ_covs %*% t(grab(mc, "beta_occ")) +
		  my_data$cell_area
		
		# convert it to a probability
		latent <- 1 - exp(-exp(latent))
		
		# Now calculate the conflict potential
		thin_prob <- 
			my_data$po_det_covs %*% t(grab(mc, "beta_po_det")) +
			matrix(
				grab(mc,"po_mu") + 
			  grab(mc, paste0("po_season\\[",year,"\\]")),
				ncol = nrow(mc),
				nrow = my_data$G,
				byrow = TRUE
		)
		
		# convert it to a probability
		thin_prob <- plogis(thin_prob)
		
		# take their product to get the conditional conflict estimate
		#  for the presence only points test data.
		this_yrs_conflict <- my_data$po_pixel[my_data$opp_year == year]
		con_conflict[[year]] <- 
			latent[this_yrs_conflict,] *
			thin_prob[this_yrs_conflict,]
	}
	# combine all the data together now
	all_con <- do.call("rbind", con_conflict)
	
	# calculate the brier score for each of these data points and posterior
	#  simulations
	
	brier[,sp] <- (1 / nrow(all_con)) * (colSums((all_con - 1)^2))
}

apply(brier, 2, quantile, probs = c(0.025,0.5,0.975))

saveRDS(
	brier,
	"./mcmc_output/validation/brier_score.RDS"
)

set.seed(444)
my_samp <- sample(1:25000, 5000)


# Calculate ROC and AUC
sp_list <- vector("list", length = 3)
for(sp in 1:3){
	species <- names(my_sum)[sp]
	source("./R/format_data_for_validation.R")
	
	# generate predictions for the validation data
	mc <- do.call("rbind", my_mcmc[[sp]]$mcmc)
	mc <- mc[my_samp,]
	thresh_list<- vector("list", length = my_data$nyear)
	pb <- txtProgressBar(max = my_data$nyear)
	for(year in 1:my_data$nyear){
		setTxtProgressBar(pb, year)
		# Get latent occupancy across the city
		latent <- my_data$X %*% t(grab(mc, paste0("b\\[\\w\\w?,",year,"\\]"))) +
			my_data$occ_covs %*% t(grab(mc, "beta_occ")) +
			my_data$cell_area
		
		# convert it to a probability
		latent <- 1 - exp(-exp(latent))
		
		# Now calculate the conflict potential
		thin_prob <- 
			my_data$po_det_covs %*% t(grab(mc, "beta_po_det")) +
			matrix(
				grab(mc,"po_mu") + 
					grab(mc, paste0("po_season\\[",year,"\\]")),
				ncol = nrow(mc),
				nrow = my_data$G,
				byrow = TRUE
			)
		
		# convert it to a probability
		thin_prob <- plogis(thin_prob)
		
		# take their product to get the conditional conflict estimate
		#  for the presence only points test data.
		my_est <- latent * thin_prob
		
		this_yrs_conflict <- my_data$po_pixel[my_data$opp_year == year]
		#con_conflict[[year]] <- 
		#	my_est
		
		my_thresh <- seq(0,1,0.025)
		TPRmat <- FPRmat <- ACCmat <- matrix(
			NA,
			nrow = length(my_thresh),
			ncol = ncol(my_est)
		)

		
		for(i in 1:length(my_thresh)){
			pred <- my_est > my_thresh[i]
			TP <- FP <-  pred == 1 
			TP <- colSums(TP[this_yrs_conflict,])
			FP <- colSums(FP[-this_yrs_conflict,])
			TN <- FN <- pred == 0
			TN <- colSums(TN[-this_yrs_conflict,])
			FN <- colSums(FN[this_yrs_conflict,])
			TPRmat[i,] <- TP / (TP + FN)
			FPRmat[i,] <- FP / (FP + TN)
			ACCmat[i,] <- (TP + TN) / (TP + TN + FP + FN)
		}
		thresh_list[[year]] <- list(
			tpr = TPRmat,
			fpr = FPRmat,
			acc = ACCmat
		)
		
	}
	sp_list[[sp]] <- thresh_list

}
names(sp_list) <- names(my_sum)
saveRDS(
	sp_list,
	"./mcmc_output/validation/ROC_scores.RDS"
)

# from here calculate AUC (overall, plus for each season).

# calculates the AUC
area<-function(x,y){
	X <- cbind(x,y)
	X<-rbind(X,X[1,])
	x<-X[,1]
	y<-X[,2] 
	lx<-length(x)
	abs(sum((x[2:lx]-x[1:lx-1])*(y[2:lx]+y[1:lx-1]))/2)
}
auc_list <- vector("list", length = 3)
for(sp in 1:3){
	
	# get mcmc together for a given species
	fpr <- tpr <- acc <-  array(
		NA,
		dim = c(
			length(my_thresh),
			length(my_samp),
			my_data$nyear
		)
	)
	
	for(year in 1:my_data$nyear){
		tpr[,,year] <- sp_list[[sp]][[year]]$tpr
		fpr[,,year] <- sp_list[[sp]][[year]]$fpr
		acc[,,year] <- sp_list[[sp]][[year]]$acc
		
	}
		
		# ROC for each year
		tpr_med <- apply(
			tpr,
			c(1,3),
			median
		)
		
		fpr_med <- apply(
			fpr,
			c(1,3),
			median
		)
		acc_med <- apply(
			acc,
			c(1,3),
			median
		)
	
auc <- matrix(NA, nrow = dim(tpr)[2], ncol = dim(tpr)[3])
	# now go through and calculate AUC
	for(i in 1:length(my_samp)){
		for(year in 1:my_data$nyear){
			auc[i,year] <- area(
				c(0, rev(fpr[,i,year]), 1),
				c(0, rev(tpr[,i,year]), 0)
			)
		}
	}

auc_list[[sp]] <- list(
	auc_global = quantile(auc, probs = c(0.025,0.5,0.975)),
	auc_year = apply(auc, 2,quantile, probs = c(0.025,0.5,0.975)),
	tpr = tpr_med,
	fpr = fpr_med,
	acc = acc_med
)
}

names(auc_list) <- names(my_sum)

saveRDS(
	auc_list,
	"./mcmc_output/validation/model_auc.RDS"
)
