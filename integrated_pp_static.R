model{
	for(g in 1:G){
		# log-linear predictor intensity across landscape
		#  Need to calculate for each cell because we need
		#  the background estimate for presence only data
		log(lambda[g]) <- inprod(occ_covs[g,], beta_occ) + cell_area
		# logit-linear predictor for thinned Poisson process. 
		#  Similar to lambda, we need to calculate this for each
		#  cell for the background estimate.
		logit(thin_prob[g]) <- inprod(po_det_covs[g,] , beta_po_det)
		thin_lambda[g] <- lambda[g] * thin_prob[g]
		# Probability of occupancy in a cell via inverse complementary log log link
		psi[g] <- 1 - exp(-lambda[g])
		# latent occupancy in each cell. For PA data.
		z[g] ~ dbern(psi[g])
	}
	# Approximation of background interval by Riemann sum. Divide by number of PO
	#  data points so we can calculate likelihood for each PO via ones trick. 
	background <- sum(thin_lambda[]) / npo
	#b2 <- thin_prob[1]
	#b3 <- lambda[1]
	#b4 <- exp(lambda[1]) * ilogit(thin_prob[1])
	#
	# Analyze the (opp)orntunistic PO data
	for(opp in 1:npo){
  # lambda for PO
		log(lambda_po[opp]) <- inprod(occ_covs[po_pixel[opp],], beta_occ)
	# thin for PO
		logit(thin_lambda_po[opp]) <- inprod(po_det_covs[po_pixel[opp],] , beta_po_det)
	# log-likelihood for PO data
		ll_po[opp] <- log(lambda_po[opp] * thin_lambda_po[opp]) -
			            log(background)
	# exponentiate ll_po and divide by a large constant for bernoulli ones trick
	  exp_po[opp] <- exp(ll_po[opp]) / CONSTANT
	# one's trick
	  ones[opp] ~ dbern(exp_po[opp])
	}
	
	for(site in 1:npa){
		logit(pa_det_prob[site]) <- inprod(pa_det_covs[pa_pixel[site],],
													beta_pa_det)
		pa_mu[site] <- z[pa_pixel[site]] * pa_det_prob[site]
		y_pa[site] ~ dbin(pa_mu[site], 4)
	}
	
	beta_occ[1] ~ dlogis(0, 1)
	beta_occ[2] ~ dlogis(0, 1)
	beta_pa_det[1] ~ dlogis(0, 1)
	beta_po_det[1] ~ dlogis(0, 1)
	beta_po_det[2] ~ dlogis(0, 1)
	
}


