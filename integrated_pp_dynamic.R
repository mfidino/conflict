model{
	for(g in 1:G){
		for(year in 1:nyear){
		# log-linear predictor intensity across landscape
		#  Need to calculate for each cell because we need
		#  the background estimate for presence only data
		log(lambda[g, year]) <- inprod(occ_covs[g, ], beta_occ[ , year]) + cell_area
		# logit-linear predictor for thinned Poisson process. 
		#  Similar to lambda, we need to calculate this for each
		#  cell for the background estimate.
		logit(thin_prob[g, year]) <- inprod(po_det_covs[g, ] , beta_po_det[ , year])
		thin_lambda[g, year] <- lambda[g, year] * thin_prob[g, year]
		# Probability of occupancy in a cell via inverse complementary log log link
		psi[g, year] <- 1 - exp(-lambda[g, year])
		# latent occupancy in each cell. For PA data.
		z[g, year] ~ dbern(psi[g, year])
		}
	}
	# Approximation of background interval by Riemann sum. Divide by number of PO
	#  data points so we can calculate likelihood for each PO via ones trick. 
	for(year in 1:nyear){
	background[year] <- sum(thin_lambda[ , year]) / npo[year]
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
		log(lambda_po[opp]) <- inprod(occ_covs[po_pixel[opp], ], 
																	beta_occ[ , opp_year[opp]])
		# thin for PO, nested indexing
		logit(thin_lambda_po[opp]) <- inprod(po_det_covs[po_pixel[opp], ], 
																							 beta_po_det[ , opp_year[opp]])
		# log-likelihood for PO data, nested indexing
		ll_po[opp] <- log(lambda_po[opp] * thin_lambda_po[opp]) -
											log(background[opp_year[opp]])
		# exponentiate ll_po and divide by a large constant for bernoulli ones trick
		exp_po[opp] <- exp(ll_po[opp]) / CONSTANT
		# one's trick
		ones[opp] ~ dbern(exp_po[opp])
	}
	
	for(site in 1:npa){
		for(year in 1:nyear){
		logit(pa_det_prob[site, year]) <- inprod(pa_det_covs[pa_pixel[site],],
																						 beta_pa_det[ , year])
		pres_abs_mu[site, year] <- z[pa_pixel[site], year] * pa_det_prob[site, year]
		y_pa[site, year] ~ dbin(pres_abs_mu[site, year], J[site,year])
		}
	}
	# priors for latent state
	# temporal random effect
	psi_mu ~ dlogis(0,1)
	psi_tau_season ~ dgamma(1,1)
	psi_sd_season <- 1 / sqrt(psi_tau_season)
	lambda_beta_occ ~ dunif(0.001,10)
	for(year in 1:nyear){
		psi_season[year] ~ dnorm(psi_mu, psi_tau_season)
		beta_occ[1, year] <- psi_season[year] 
	}
	for(latent in 2:nlatent){
			beta_occ_fill[latent-1] ~ ddexp(0, lambda_beta_occ)
		for(year in 1:nyear){
			beta_occ[latent,year] <- beta_occ_fill[latent-1]
		}
	}
	# priors for presence absence observation 
	# temporal random effect
	pa_mu ~ dlogis(0,1)
	pa_tau_season ~ dgamma(1,1)
	pa_sd_season <- 1 / sqrt(pa_tau_season)
	lambda_pa_det ~ dunif(0.001,10)
	for(year in 1:nyear){
		pa_season[year] ~ dnorm(pa_mu, pa_tau_season)
		beta_pa_det[1, year] <- pa_season[year] 
	}
	for(obs_pa in 2:nobs_pa){
		beta_observation[obs_pa-1] ~ ddexp(0, lambda_pa_det)
		for(year in 1:nyear){
			beta_pa_det[obs_pa, year] <- beta_observation[obs_pa-1]
		}
	}
	# priors for presence only observation
	# temporal random effect
	po_mu ~ dlogis(0,1)
	po_tau_season ~ dgamma(1,1)
	po_sd_season <- 1 / sqrt(po_tau_season)
	lambda_po_det ~ dunif(0.001, 10)
	for(year in 1:nyear){
		po_season[year] ~ dnorm(po_mu, po_tau_season)
		beta_po_det[1, year] <- po_season[year]
	}
	for(obs_po in 2:nobs_po){
		beta_po_fill[obs_po-1]  ~ ddexp(0, lambda_po_det)
		for(year in 1:nyear){
			beta_po_det[obs_po, year] <- beta_po_fill[obs_po-1]
		}
	}
	
}


