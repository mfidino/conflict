#############################################
#
#
# fit the model to each species
#
# Written by M. Fidino
#
#############################################

source("sourcer.R")
packs <- c("lubridate", "raster", "sp", "sf", "runjags", "coda", "nimble")
package_load(packs)


my_species <- c("opossum", "coyote", "raccoon")

for(animal in 2:length(my_species)){

# species is used in "format_data_for_analysis.R" to 
#  pull in the correct data.
species <- my_species[animal]

# print it out to console so I can see where we are at.
cat(species)
source("format_data_for_analysis.R")

test <- readBUGSmodel("integrated_pp_dynamic_faster_ranef.R",
											data = my_data)

tmp <- my_inits(1)
longshot <- function(seed,data){
	library(nimble)
	
	my_inits <- function(chain){
		gen_list <- function(chain = chain){
			list( 
				z = matrix(1, ncol = data$nyear, nrow = data$G),
				beta_occ = rnorm(data$nlatent, 0, 0.25),
				beta_pa_det = rnorm(data$nobs_pa, 0, 0.25),
				beta_po_det = rnorm(data$nobs_po, 0, 0.25),
				psi_mu = rnorm(1, -5, 0.25),
				pa_mu = rnorm(1, -2.75, 0.25),
				po_mu = rnorm(1, 3, 0.25),
				lambda_beta_occ = runif(1, 1, 2),
				lambda_pa_det = runif(1, 1, 2),
				lambda_po_det = runif(1, 1, 2),
				psi_tau_season = rgamma(1, 2, 4),
				pa_tau_season = rgamma(1, 2, 4),
				po_tau_season = rgamma(1, 2, 4),
				psi_season = rnorm(data$nyear, 0, 0.25),
				pa_season = rnorm(data$nyear, 0, 0.25),
				po_season = rnorm(data$nyear,0, 0.25),
				theta = rnorm(1),
				.RNG.name = switch(chain,
													 "1" = "base::Wichmann-Hill",
													 "2" = "base::Wichmann-Hill",
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
	
	
	test <- readBUGSmodel("integrated_pp_dynamic_faster_ranef.R",
												data = data,
												inits = my_inits(1),
												calculate = TRUE)
	outA <- nimbleMCMC(
		code = test,
		constants = data,
		inits = my_inits(1),
		monitors = c(
			"beta_occ", "beta_pa_det", "beta_po_det",
			"lambda_beta_occ", "lambda_pa_det", "lambda_po_det",
			"psi_mu", "pa_mu", "po_mu", "psi_season",
			"pa_season", "po_season", "psi_sd_season",
			"pa_sd_season", "po_sd_season"
		),
		nburnin = 5000,
		niter = 10000,
		nchains = 1,
		progress = TRUE,
		setSeed = seed
	)
	return(outA)
	
}


library(parallel)
this_cluster <- makeCluster(3)

chain_output <- parLapply(cl = this_cluster, X = 1:3, 
													fun = longshot, 
													data = my_data)

# It's good practice to close the cluster when you're done with it.
stopCluster(this_cluster)


# Note: my_data is the data list that is created from the 
#       script above.
m1 <- run.jags(model = "integrated_pp_dynamic_faster_ranef.R", 
							 data = my_data, 
							 n.chains = 3, 
							 inits = inits, 
							 monitor = c(
							 	"beta_occ", "beta_pa_det", "beta_po_det",
							 	"lambda_beta_occ", "lambda_pa_det", "lambda_po_det",
							 	"psi_mu", "pa_mu", "po_mu", "psi_season",
							 	"pa_season", "po_season", "psi_sd_season",
							 	"pa_sd_season", "po_sd_season"
							 ), 
							 adapt = 100, 
							 burnin = 2000, 
							 sample = 1000, 
							 thin = 1,
							 modules = "glm",
							 method = 'parallel',
							 summarise = FALSE)

# check for convergence and the like

msum <- round(summary(m1),2)
saveRDS(
	m1,
	paste0("./mcmc_output/",species,"_model_1mo.RDS")
)

# save plots of the chains to inspect them
mc <- as.matrix(as.mcmc.list(m1), chains = TRUE)

to_write <- colnames(mc)
for(i in 2:ncol(mc)){
	jpeg(
		paste0(
			"./mcmc_output/diagnostic_plots/",
			species, "/",to_write[i],".jpeg")
	)
	plot(mc[,i] ~ c(1:nrow(mc)), type = 'p',
			 main = to_write[i], col = rainbow(5)[mc[,1]])	
	dev.off()
}

}


# run raccoon model longer
m1 <- readRDS("./mcmc_output/raccoon_model.RDS")

m1 <- extend.jags(m1, adapt = 100, sample = 10000)
