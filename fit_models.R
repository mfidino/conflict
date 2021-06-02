#############################################
#
#
# fit the model to each species
#
# Written by M. Fidino
#
#############################################




my_species <- c("raccoon", "opossum", "coyote")

for(animal in 1:length(my_species)){
	my_species <- c("opossum", "raccoon", "coyote")
	source("sourcer.R")
	packs <- c("lubridate", "raster", "sp", "sf", "runjags", "coda", "mgcv")
	package_load(packs)

# species is used in "format_data_for_analysis.R" to 
#  pull in the correct data.
species <- my_species[animal]

# print it out to console so I can see where we are at.
cat(species)
source("format_data_for_analysis.R")

#y_pa_init <- my_data$y_pa
#y_pa_init[is.na(y_pa_init)] <- 0
#y_pa_init[!is.na(my_data$y_pa)] <- NA

# my_inits <- function(chain){
# 	gen_list <- function(chain = chain){
# 		list(
# 			z = matrix(1, ncol = my_data$nyear, nrow = my_data$G),
# 			beta_occ = rnorm(my_data$nlatent, 0, 0.25),
# 			temp_occ = rnorm(2),
# 			beta_pa_det = rnorm(my_data$nobs_pa, 0, 0.25),
# 			beta_po_det = rnorm(my_data$nobs_po, 0, 0.25),
# 			beta_po_det_mu = rnorm(my_data$nobs_po),
# 			beta_po_det_tau = rgamma(my_data$nobs_po,1,1),
# 			psi_mu = rnorm(1, -5, 0.25),
# 			pa_mu = rnorm(1, -2.75, 0.25),
# 			po_mu = rnorm(1, 3, 0.25),
# 			g_tau = rgamma(1, 2, 4),
# 			g_re = rnorm(my_data$G, 0, 0.25),
# 			pa_tau_season = rgamma(1, 2, 4),
# 			po_tau_season = rgamma(1, 2, 4),
# 			pa_season = rnorm(my_data$nyear, 0, 0.25),
# 			po_season = rnorm(my_data$nyear,0, 0.25),
# 			y_pa = y_pa_init,
# 			.RNG.name = switch(chain,
# 												 "1" = "base::Wichmann-Hill",
# 												 "2" = "base::Wichmann-Hill",
# 												 "3" = "base::Super-Duper",
# 												 "4" = "base::Mersenne-Twister",
# 												 "5" = "base::Wichmann-Hill",
# 												 "6" = "base::Marsaglia-Multicarry",
# 												 "7" = "base::Super-Duper",
# 												 "8" = "base::Mersenne-Twister"),
# 			.RNG.seed = sample(1:1e+06, 1)
# 		)
# 	}
# 	return(switch(chain,
# 								"1" = gen_list(chain),
# 								"2" = gen_list(chain),
# 								"3" = gen_list(chain),
# 								"4" = gen_list(chain),
# 								"5" = gen_list(chain),
# 								"6" = gen_list(chain),
# 								"7" = gen_list(chain),
# 								"8" = gen_list(chain)
# 	)
# 	)
# }

# Note: my_data is the data list that is created from the 
#       script above.
m1 <- run.jags(model = "integrated_pp_dynamic_prd_sre.R", 
							 data = my_data, 
							 n.chains = 3, 
							 inits = my_inits, 
							 monitor = c(
							 	"beta_occ", "beta_pa_det", "beta_po_det",
							 	"psi_mu", "pa_mu", "po_mu",
							 	"pa_season", "po_season", "b", "rho",
							 	"gam_sd",
							 	"pa_sd_season", "po_sd_season"
							 ), 
							 adapt = 1000, 
							 burnin = 10000, 
							 sample = 2000, 
							 thin = 4,
							 modules = "glm",
							 method = 'parallel',
							 summarise = FALSE)

# check for convergence and the like

msum <- round(summary(m1),2)
saveRDS(
	m1,
	paste0("./mcmc_output/",species,"_model.RDS")
)

#saveRDS(
#	m2,
#	paste0("./mcmc_output/coyote_theta_model.RDS")
#)



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
rm(list = ls())

}


# run raccoon model longer
m1 <- readRDS("./mcmc_output/coyote_model.RDS")

m3 <- extend.jags(m2, adapt = 100, sample = 3000)
