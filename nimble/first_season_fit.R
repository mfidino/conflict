one_season <- function(mseed, my_data){
	library(nimble)
	y_pa_init <- my_data$y_pa
	y_pa_init[is.na(y_pa_init)] <- 0
	y_pa_init[!is.na(my_data$y_pa)] <- NA
	
	my_inits <- function(chain){
		gen_list <- function(chain = chain){
			list( 
				z = matrix(1, ncol = my_data$nyear, nrow = my_data$G),
				beta_occ = rnorm(my_data$nlatent, 0, 0.25),
				beta_pa_det = rnorm(my_data$nobs_pa, 0, 0.25),
				beta_po_det = rnorm(my_data$nobs_po, 0, 0.25),
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
	the_mod <- nimbleCode({
		# I'm writing this model in a way that is slightly less comprehensible
		#  But it runs about 20% faster.
		for(g in 1:G){
			for(year in 1:nyear){
				# log-linear predictor intensity across landscape
				#  Need to calculate for each cell because we need
				#  the background estimate for presence only data
				log(lambda[g, year]) <- psi_mu + cell_area + 
					inprod(occ_covs[g,1:nlatent ], beta_occ[1:nlatent]) 
				# logit-linear predictor for thinned Poisson process. 
				#  Similar to lambda, we need to calculate this for each
				#  cell for the background estimate.
				logit(thin_prob[g, year]) <- po_mu +
					inprod(po_det_covs[g,1:nlatent ] , beta_po_det[1:nlatent])
				# Probability of occupancy in a cell via inverse complementary log log link
				# latent occupancy in each cell. For PA data.
				z[g, year] ~ dbern(1 - exp(-lambda[g, year]))
			}
		}
		# Approximation of background interval by Riemann sum. Divide by number of PO
		#  data points so we can calculate likelihood for each PO via ones trick. 
		for(year in 1:nyear){
			background[year] <- inprod(lambda[1:G, year], thin_prob[1:G, year] ) / npo
		}
		#
		# Analyze the (opp)orntunistic PO data
		#  Since there is varying amounts of data per year
		#  we use nested indexing to collect the appropriate parameters
		#  for each PO datapoint.
		#  This means that we have grouped all of the PO data into one long vector.
		#  of length(all_npo)
		for(opp in 1:all_npo){
			# lambda for PO, nested indexing
			log(lambda_po[opp]) <- psi_mu + cell_area + 
				inprod(occ_covs[po_pixel[opp],1:nlatent ], beta_occ[1:nlatent])
			# thin for PO, nested indexing
			logit(thin_lambda_po[opp]) <- po_mu +  
				inprod(po_det_covs[po_pixel[opp], 1:nobs_po], beta_po_det[1:nobs_po])
			# one's trick
			ones[opp] ~ dbern(
				exp(
					log(lambda_po[opp] * thin_lambda_po[opp]) -
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
		psi_mu ~ T(dnorm(-3.2, 1), -8.5, 5)
		for(latent in 1:nlatent){
			beta_occ[latent] ~  dnorm(0, 0.45)
		}
		# priors for presence absence observation 
		# temporal random effect
		pa_mu ~ T(dnorm(0, 0.75), -5,5)
		for(obs_pa in 1:nobs_pa){
			beta_pa_det[obs_pa] ~ dnorm(0, 0.75)
		}
		# priors for presence only observation
		# temporal random effect
		po_mu ~ T(dnorm(0, 0.75),-5,5)
		for(obs_po in 1:nobs_po){
			beta_po_det[obs_po]  ~  dnorm(0, 0.75)
		}
	})
	set.seed(seed = mseed )
	tc <- my_data[
		c("G", "pa_pixel", "po_pixel",
			"opp_year", "npo", "all_npo", "cell_area", "CONSTANT", "nlatent",
			"nobs_po", "nobs_pa", "nyear", "npa")]
	td <- my_data[-which(names(my_data) %in% c("G", "pa_pixel", "po_pixel",
																						 "opp_year", "npo", "all_npo", "cell_area", "CONSTANT", "nlatent",
																						 "nobs_po", "nobs_pa", "nyear", "npa"))]
	
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
			"psi_mu", "pa_mu", "po_mu", "lambda[1,1]"
		),
		control = list(maxContractionsWarning = FALSE)
	)
	conf$addSampler(
		target = c(
			"beta_occ", "beta_pa_det", "beta_po_det",
			"psi_mu", "pa_mu", "po_mu"
		),
		type = "RW_block",
		control = list(adaptInterval = 10000)
	)
		

	Rmcmc <- nimble::buildMCMC(conf = conf)
	
	
	Cmodel <- compileNimble(the_model)
	Cmcmc <- compileNimble(
		Rmcmc,
		project = the_model,
		resetFunctions = TRUE
	)
	
	Cmcmc$run(10, thin = 1)
	
	samples2 <- as.matrix(Cmcmc$mvSamples)
	
	tc <- samples2[,13]
	ay <- samples2[1,c(ncol(samples2),1:5)]
	aa <- my_data$occ_covs[1,]
	
	tt <- ay %*% c(1, aa)
	return(samples)
	
}

# here is a function to visually inspect output, and caluculate
#  the gelman rubin diagnostic
inspect_fit <- function(x){
	
}

