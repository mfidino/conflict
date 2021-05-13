#####################################
#
# Adding a temporal term to the model
#
#####################################

my_species <- c("coyote", "opossum", "raccoon")

source("sourcer.R")

packs <- c("lubridate", "raster", "sp", "sf", "runjags", "coda", "nimble")
package_load(packs)

for(animal in 1:length(my_species)){
	species <- my_species[animal]
	source("format_data_for_analysis.R")
	

	
	# specify the temporal model
	the_mod <- nimbleCode({
		# I'm writing this model in a way that is slightly less comprehensible
		#  But it runs about 20% faster.
		for(g in 1:G){
			for(year in 1:nyear){
				# log-linear predictor intensity across landscape
				#  Need to calculate for each cell because we need
				#  the background estimate for presence only data
				log(lambda[g, year]) <- psi_mu + cell_area + 
					inprod(occ_covs[g,1:nlatent ], beta_occ[1:nlatent]) #+
				#	inprod(temp_covs[year,1:2], temp_occ[1:2])
				# logit-linear predictor for thinned Poisson process. 
				#  Similar to lambda, we need to calculate this for each
				#  cell for the background estimate.
				logit(thin_prob[g, year]) <- po_mu + po_year[year] +
					inprod(po_det_covs[g,1:nlatent ] , beta_po_det[1:nlatent])
				# Probability of occupancy in a cell via inverse complementary log log link
				# latent occupancy in each cell. For PA data.
				z[g, year] ~ dbern(1 - exp(-lambda[g, year]))
			}
		}
		# Approximation of background interval by Riemann sum. Divide by number of PO
		#  data points so we can calculate likelihood for each PO via ones trick. 
		for(year in 1:nyear){
			background[year] <- inprod(lambda[1:G, year], thin_prob[1:G, year] ) / npo[year]
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
				logit(pa_det_prob[site, year]) <- pa_mu + pa_year[year]+
					inprod(pa_det_covs[pa_pixel[site],1:nobs_pa],beta_pa_det[1:nobs_pa])
				y_pa[site, year] ~ dbin(
					z[pa_pixel[site], year] * pa_det_prob[site, year],
					J[site,year]
				)
			}
		}
		# priors for latent state
		# temporal random effect
		psi_mu ~ dnorm(-3.2, 1.3)
		for(latent in 1:nlatent){
			beta_occ[latent] ~  dnorm(0, 1.3)
		}
		#for(btime in 1:2){
		#	temp_occ[btime] ~ dnorm(0, 1.3)
		#}
		#shape ~ dunif(0.001, 20)
		#scale ~ dunif(0.001, 20)
		#for(year in 1:nyear){
	  #		amp[year] ~ dgamma(shape, scale)
		#}
		# priors for presence absence observation 
		# temporal random effect
		pa_mu ~ dnorm(0, 1.3)
		for(obs_pa in 1:nobs_pa){
			beta_pa_det[obs_pa] ~ dnorm(0, 1.3)
		}
		# temporal random effects for detection level
		#   of the model.
		pa_tau ~ dgamma(1,1)
		pa_sd <- 1 / sqrt(pa_tau)
		po_tau ~ dgamma(1,1)
		po_sd <- 1 / sqrt(po_tau)
		for(year in 1:nyear){
			pa_year[year] ~ dnorm(0, pa_tau)
			po_year[year] ~ dnorm(0, po_tau)
		}
		# priors for presence only observation
		# temporal random effect
		po_mu ~ dnorm(0, 1.3)
		for(obs_po in 1:nobs_po){
			beta_po_det[obs_po]  ~  dnorm(0, 1.3)
		}
	})
	
	
	the_mod <- nimbleCode({
		# I'm writing this model in a way that is slightly less comprehensible
		#  But it runs about 20% faster.
		for(g in 1:G){
			for(year in 1:nyear){
				y1[g, year] ~ dbern(psi1[g,year])
				cloglog(psi1[g,year]) <- cell_area + log(lambda1[g,year]) +
					log(thin_prob[g,year])
				log(lambda1[g,year]) <- psi_mu + 
					inprod(occ_covs[g,1:nlatent], beta_occ[1:nlatent])
				logit(thin_prob[g, year]) <- po_mu + po_year[year] +
					inprod(po_det_covs[g,1:nlatent ] , beta_po_det[1:nlatent])
			}
		}
		for(site in 1:npa){
			for(year in 1:nyear){
				z[site,year] ~ dbern(psi2[site,year])
				cloglog(psi2[site,year]) <- cell_area + log(lambda2[site,year])
				log(lambda2[site,year]) <- psi_mu +
					inprod(
						occ_covs[pa_pixel[site], 1:nlatent],
						beta_occ[1:nlatent]
					)
				logit(pa_det_prob[site, year]) <- pa_mu + pa_year[year]+
					inprod(pa_det_covs[pa_pixel[site],1:nobs_pa],beta_pa_det[1:nobs_pa])
				y_pa[site, year] ~ dbin(
					z[site, year] * pa_det_prob[site, year],
					J[site,year]
				)
			}
		}
		# priors for latent state
		# temporal random effect
		psi_mu ~ dnorm(-3.2, 1.3)
		for(latent in 1:nlatent){
			beta_occ[latent] ~  dnorm(0, 1.3)
		}
		#for(btime in 1:2){
		#	temp_occ[btime] ~ dnorm(0, 1.3)
		#}
		#shape ~ dunif(0.001, 20)
		#scale ~ dunif(0.001, 20)
		#for(year in 1:nyear){
		#		amp[year] ~ dgamma(shape, scale)
		#}
		# priors for presence absence observation 
		# temporal random effect
		pa_mu ~ dnorm(0, 1.3)
		for(obs_pa in 1:nobs_pa){
			beta_pa_det[obs_pa] ~ dnorm(0, 1.3)
		}
		# temporal random effects for detection level
		#   of the model.
		pa_tau ~ dgamma(1,1)
		pa_sd <- 1 / sqrt(pa_tau)
		po_tau ~ dgamma(1,1)
		po_sd <- 1 / sqrt(po_tau)
		for(year in 1:nyear){
			pa_year[year] ~ dnorm(0, pa_tau)
			po_year[year] ~ dnorm(0, po_tau)
		}
		# priors for presence only observation
		# temporal random effect
		po_mu ~ dnorm(0, 1.3)
		for(obs_po in 1:nobs_po){
			beta_po_det[obs_po]  ~  dnorm(0, 1.3)
		}
	})
	
	ack <- my_data
	# make a matrix for y1
	y1 <- matrix(0, ncol = ack$nyear, nrow = ack$G)
	for(i in  1:ncol(y1)){
		to_add <- which(ack$opp_year == i)
		y1[ack$po_pixel[to_add],i] <- 1
	}
	ack$y1 <- y1
	
	tc <- my_data[
		c("G", "pa_pixel", "po_pixel",
			"opp_year", "npo", "all_npo", "cell_area", "CONSTANT", "nlatent",
			"nobs_po", "nobs_pa", "nyear", "npa")]
	
	td <- my_data[-which(names(my_data) %in% c("G", "pa_pixel", "po_pixel",
																						 "opp_year", "npo", "all_npo", "cell_area", "CONSTANT", "nlatent",
																						 "nobs_po", "nobs_pa", "nyear", "npa"))]
	td$y1 <- y1
	
	# generate temporal terms
	cov1 <- sin(2 * pi * (1:my_data$nyear )/4)
	cov2 <- cos(2 * pi * (1:my_data$nyear )/4)
	td$temp_covs <- matrix(c(cov1, cov2), ncol = 2, nrow = length(cov1))
	
	my_inits <- function(chain){
		gen_list <- function(chain = chain){
			list(
				z = matrix(1, ncol = my_data$nyear, nrow = my_data$npa),
				beta_occ = rnorm(my_data$nlatent, 0, 0.25),
				beta_pa_det = rnorm(my_data$nobs_pa, 0, 0.25),
				beta_po_det = rnorm(my_data$nobs_po, 0, 0.25),
				temp_occ = rnorm(2),
				shape = runif(1, 0.5, 1.5),
				scale = runif(1, 0.5, 1.5),
				amp = runif(my_data$nyear, 0.9, 1.1),
				psi_mu = rnorm(1, -5, 0.25),
				pa_mu = rnorm(1, -2.75, 0.25),
				po_mu = rnorm(1, 3, 0.25),
				pa_tau = rgamma(1, 2, 4),
				po_tau = rgamma(1, 2, 4),
				pa_year = rnorm(my_data$nyear, 0, 0.25),
				po_year = rnorm(my_data$nyear,0, 0.25),
				y_pa = y_pa_init)
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
	
	
	the_model <- nimbleModel(
		the_mod,
		constants = tc,
		data = td,
		inits = my_inits(1)
	)
	conf <- configureMCMC(
		the_model,
		monitors = c(
			"beta_occ", "beta_pa_det", "beta_po_det",
			"psi_mu", "pa_mu", "po_mu",
			"po_sd", "pa_sd", "po_year", "pa_year"
		),
		control = list(maxContractionsWarning = FALSE)
	)
	conf$addSampler(
		target = c(
			"beta_occ", "beta_pa_det", "beta_po_det",
			"psi_mu", "pa_mu", "po_mu",
		  "po_year", "pa_year"
		),
		type = "RW_block",
		control = list(maxContractionsWarning = FALSE)
	)
	
	
	Rmcmc <- nimble::buildMCMC(conf = conf)
	
	
	Cmodel <- compileNimble(the_model)
	Cmcmc <- compileNimble(
		Rmcmc,
		project = the_model,
		resetFunctions = TRUE
	)
	
	Cmcmc$run(5000, reset = FALSE)
	
	samples <- as.data.frame(as.matrix(Cmcmc$mvSamples))
	
	plot(samples$pa_mu, type = 'l')
	plot(samples$`temp_occ[2]`, type = 'l')
	
	}