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
	#  Since there is varying amounts of data per years
	#  we use nested indexing to collect the appropriate parameters
	#  for each PO datapoint.
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
		pa_mu[site, year] <- z[pa_pixel[site], year] * pa_det_prob[site, year]
		y_pa[site, year] ~ dbin(pa_mu[site, year], 4)
		}
	}
	# priors for latent state
	for(latent in 1:nlatent){
		for(year in 1:nyear){
			beta_occ[latent, year] ~ dlogis(0, 1)
		}
	}
	# priors for presence absence observation 
	#  currently does not vary by year
	for(obs_pa in 1:nobs_pa){
		beta_observation[obs_pa] ~ dlogis(0, 1)
		for(year in 1:nyear){
			beta_pa_det[obs_pa, year] <- beta_observation[obs_pa]
		}
	}
	# priors for presence only observation
	#  currently does not vary by year
	for(obs_po in 1:nobs_po){
		beta_po_fill[obs_po]  ~ dlogis(0, 1)
		for(year in 1:nyear){
			beta_po_det[ops_po, year] <- beta_po_fill[obs_po]
		}
	}
	
}


