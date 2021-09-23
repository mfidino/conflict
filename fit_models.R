#############################################
#
#
# fit the model to each species
#
# Written by M. Fidino
#
#############################################

my_species <- c("opossum", "raccoon", "coyote")

for(animal in 1:length(my_species)){
	
	source("sourcer.R")
	packs <- c(
		"lubridate", "raster", "sp", "sf", "runjags", "coda", "mgcv"
	)
	package_load(packs)

  # species is used in "format_data_for_analysis.R" to 
  #  pull in the correct data.
  species <- my_species[animal]
  
  # print it out to console so I can see where we are at.
  cat(species)
  source("format_data_for_analysis.R")
  
  # Note: my_data is the data list that is created from the 
  #       script above.
  m1 <- run.jags(model = "./JAGS/dynamic_integrated_occupancy_gam.R", 
  							 data = my_data, 
  							 n.chains = 4, 
  							 inits = my_inits, 
  							 monitor = c(
  							 	"beta_occ", "beta_pa_det", "beta_po_det",
  							 	"psi_mu", "pa_mu", "po_mu",
  							 	"pa_season", "po_season", "b", "rho",
  							 	"gam_sd",
  							 	"pa_sd_season", "po_sd_season"
  							 ), 
  							 adapt = 1000, 
  							 burnin = 100000, 
  							 sample = 6250, 
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
  # remove everything but the my_species object
  to_go <- ls()
  to_go <- to_go[-grep("my_species", to_go)]
  rm(list = to_go)

}
