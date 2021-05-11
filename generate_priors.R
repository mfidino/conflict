#############################################
#
#
# fit the model to each species
#
# Written by M. Fidino
#
#############################################





	


# make a folder to store the prior results stuff
if(!file.exists("prior_fit")){
	dir.create("./prior_fit")
	dir.create("./prior_fit/coyote")
	dir.create("./prior_fit/opossum")
	dir.create("./prior_fit/raccoon")
}

# loop through each species, and fit the model.
for(animal in 1:3){
	my_species <- c("opossum", "coyote", "raccoon")
	
	source("sourcer.R")
	packs <- c(
		"lubridate", "raster", "sp", "runjags",
		"sf", "coda", "nimble", "parallel"
	)
	package_load(packs)
	
	y_pa_init <- first_fit$y_pa
	y_pa_init[is.na(y_pa_init)] <- 0
	y_pa_init[!is.na(first_fit$y_pa)] <- NA

	
	first_inits <- function(chain){
		gen_list <- function(chain = chain){
			list( 
				z = matrix(1, ncol = first_fit$nyear, nrow = first_fit$G),
				beta_occ = rnorm(first_fit$nlatent, 0, 0.25),
				beta_pa_det = rnorm(first_fit$nobs_pa, 0, 0.25),
				beta_po_det = rnorm(first_fit$nobs_po, 0, 0.25),
				psi_mu = rnorm(1, -5, 0.25),
				pa_mu = rnorm(1, -2.75, 0.25),
				po_mu = rnorm(1, 3, 0.25),
				y_pa = y_pa_init
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
  # species is used in "format_data_for_analysis.R" to 
  #  pull in the correct data.
  species <- my_species[animal]
  
  # print it out to console so I can see where we are at.
  cat(species)
  # and then prepare the data for analysis. Object we should know about
  # first_fit. This is the list to fit the model and generate some
  #  priors.
  source("format_data_for_analysis.R")
  
  
  
  
  m1 <- run.jags(model = "integrated_pp_static_faster_ranef.R", 
  							 data = first_fit, 
  							 n.chains = 3, 
  							 inits = first_inits, 
  							 monitor = c(
  							 	"beta_occ", "beta_pa_det", "beta_po_det",
  							 	"psi_mu", "pa_mu", "po_mu", "zsum"
  							 ), 
  							 adapt = 100, 
  							 burnin = 10000, 
  							 sample = 10000, 
  							 thin = 1,
  							 modules = "glm",
  							 method = 'parallel',
  							 summarise = FALSE)
  
  # check for convergence and the like
  
  msum <- round(summary(m1),2)
  saveRDS(
  	m1,
  	paste0("./prior_fit/",species,"_model.RDS")
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
  			"./prior_fit/",
  			species, "/",to_write[i],".jpeg")
  	)
  	plot(mc[,i] ~ c(1:nrow(mc)), type = 'p',
  			 main = to_write[i], col = rainbow(5)[mc[,1]])	
  	dev.off()
  }
	
	rm(list = ls())
}



	
	# fit in nimble
	
	
	for(animal in 2:3){
		
		my_species <- c("opossum", "coyote", "raccoon")
		
		source("sourcer.R")
		packs <- c(
			"lubridate", "raster", "sp",
			"sf", "coda", "nimble", "parallel"
		)
		package_load(packs)
		species <- my_species[animal]
		
		source("format_data_for_analysis.R")
		fit_nimble <- function(mseed, my_data){
		library(nimble)
		y_pa_init <- my_data$y_pa
		y_pa_init[is.na(y_pa_init)] <- 0
		y_pa_init[!is.na(my_data$y_pa)] <- NA
		
		first_inits <- function(chain){
			gen_list <- function(chain = chain){
				list( 
					z = matrix(1, ncol = my_data$nyear, nrow = my_data$G),
					beta = matrix(
						rnorm((my_data$nlatent + 1)*2, 0, 0.25),
						ncol = my_data$nlatent + 1,
						nrow = 2
					),
					beta_pa_det = rnorm(my_data$nobs_pa, 0, 0.25),
					pa_mu = rnorm(1, -2.75, 0.25),
					p = runif(my_data$nlatent + 1, -1, 1),
					y_pa = y_pa_init
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
	
	
	the_mod <- nimbleCode({
		# I'm writing this model in a way that is slightly less comprehensible
		#  But it runs about 20% faster.
		for(g in 1:G){
			for(year in 1:nyear){
				# log-linear predictor intensity across landscape
				#  Need to calculate for each cell because we need
				#  the background estimate for presence only data
				log(lambda[g, year]) <- cell_area + 
					inprod(occ_covs[g,1:nlatent ], beta[1,1:nlatent]) 
				# logit-linear predictor for thinned Poisson process. 
				#  Similar to lambda, we need to calculate this for each
				#  cell for the background estimate.
				logit(thin_prob[g, year]) <- 
					inprod(po_det_covs[g,1:nlatent ] , beta[2,1:nlatent])
				# Probability of occupancy in a cell via inverse complementary log log link
				# latent occupancy in each cell. For PA data.
				z[g, year] ~ dbern(1 - exp(-lambda[g, year]))
			}
		}
		# Approximation of background interval by Riemann sum. Divide by number of PO
		#  data points so we can calculate likelihood for each PO via ones trick. 
		for(year in 1:nyear){
			background[year] <- inprod(lambda[1:G, year], thin_prob[1:G, year] ) / npo
			zsum[year] <- sum(z[1:G, year])
		}
		#
		# Analyze the (opp)orntunistic PO data
		#  Since there is varying amounts of data per year
		#  we use nested indexing to collect the appropriate parameters
		#  for each PO datapoint.
		#  This means that we have grouped all of the PO data into one long vector.
		#  of length(all_npo)
		for(opp in 1:all_npo){
			# one's trick
			ones[opp] ~ dbern(
				exp(
					log(lambda[po_pixel[opp], opp_year[opp]] *
								thin_prob[po_pixel[opp], opp_year[opp]]) -
						log(background[opp_year[opp]])
				) / CONSTANT
			)
		}
		
		for(site in 1:npa){
			for(year in 1:nyear){
				logit(pa_det_prob[site, year]) <- pa_mu + 
					inprod(pa_det_covs[pa_pixel[site],1:nobs_pa],beta_pa_det[1:nobs_pa])
				y_pa[site, year] ~ dbin(
					z[pa_pixel[site], year] * pa_det_prob[site, year],
					J[site,year]
				)
			}
		}
		# priors for latent state
		# temporal random effect
		for(latent in 1:nlatent){
			p[latent] ~ dunif(-1,1)
			Sigma[1,1,latent] <- 1
			Sigma[2,2,latent] <- 1
			Sigma[2,1,latent] <- 1 * p[latent] 
			Sigma[1,2,latent] <- 1 * p[latent] 
			beta[1:2, latent] ~  dmnorm(my_means[1:2,latent], Sigma[1:2,1:2,latent])
		}
		# priors for presence absence observation 
		# temporal random effect
		pa_mu ~ T(dnorm(0, 0.75), -5,5)
		for(obs_pa in 1:nobs_pa){
			beta_pa_det[obs_pa] ~ dnorm(0, 0.75)
		}
		# priors for presence only observation
		# temporal random effect
	})
	set.seed(seed = mseed )
	tc <- my_data[
		c("G", "pa_pixel", "po_pixel",
			"opp_year", "npo", "all_npo", "cell_area", "CONSTANT", "nlatent",
			"nobs_po", "nobs_pa", "nyear", "npa")]
	td <- my_data[-which(names(my_data) %in% c("G", "pa_pixel", "po_pixel",
																						 "opp_year", "npo", "all_npo", "cell_area", "CONSTANT", "nlatent",
																						 "nobs_po", "nobs_pa", "nyear", "npa"))]
	td$my_means <- matrix(
		0,
		ncol = my_data$nlatent + 1,
		nrow = 2
	)
	tc$nlatent <- 6
	td$occ_covs <- cbind(1, td$occ_covs)
	td$po_det_covs <- cbind(1, td$po_det_covs)
	
	the_model <- nimbleModel(
		the_mod,
		constants = tc,
		data = td,
		inits = first_inits(1)
	)

	conf <- configureMCMC(
		the_model,
		monitors = c(
			"beta", "beta_pa_det",
			 "pa_mu", "zsum", "p"
		),
		control = list(maxContractionsWarning = FALSE)
	)
	for(i in 1:6){
		for(j in 1:2){
			conf$addSampler(
				target = paste0("beta[",j,",",i,"]"),
				type = "ess",
				control = (list(maxContractionsWarning = FALSE))
			)
			
		}
	}

	conf$addSampler(
		target = c( "beta_pa_det", "pa_mu"),
		type = "RW_block",
		control = list(maxContractionsWarning = FALSE)
	)
	conf$addSampler(
		target = "p",
		type = "RW_block"
	)
	
	
	Rmcmc <- nimble::buildMCMC(conf = conf)
	
	
	Cmodel <- compileNimble(the_model)
	Cmcmc <- compileNimble(
		Rmcmc,
		project = the_model,
		resetFunctions = TRUE
	)
	
	Cmcmc$run(20000)
	
	samples <- as.matrix(Cmcmc$mvSamples)
	
	return(samples)
	
		}
		
		library(parallel)
		this_cluster <- makeCluster(3)

		chain_output <- parLapply(cl = this_cluster, X = 1:3,
															fun = fit_nimble,
															my_data = first_fit)
		saveRDS(
			chain_output, 
			paste0(
				"./prior_fit/",species,"_mcmc.RDS"
			)
		)

		# It's good practice to close the cluster when you're done with it.
		stopCluster(this_cluster)
tp <- ncol(chain_output[[1]])
to_write <- colnames(chain_output[[1]])
		for(i in 1:tp){
			#dat <- sapply(
		#		chain_output,
			#	function(x)
			#		x[,i]
			#)
			dat <- chain_output
				jpeg(
					paste0(
						"./prior_fit/",
						species, "/",to_write[i],".jpeg")
				)
				plot(dat[,i] ~ c(1:nrow(dat)), type = "l",
						 ylim = range(dat[,i]), main = to_write[i])
				#lines(dat[,2] ~ c(1:nrow(dat)), col = "green")
				#lines(dat[,3] ~ c(1:nrow(dat)), col = "blue")	
				dev.off()

		}
rm(list = ls())
	}
m1 <- readRDS("./prior_fit/opossum_mcmc.RDS")


test <- as.mcmc.list(lapply(m1, as.mcmc))
		
ll <- gelman.diag(test)
