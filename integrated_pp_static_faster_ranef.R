model{
	
	# I'm writing this model in a way that is slightly less comprehensible
	#  But it runs about 20% faster.
	for(g in 1:G){
		for(year in 1:nyear){
			# log-linear predictor intensity across landscape
			#  Need to calculate for each cell because we need
			#  the background estimate for presence only data
			log(lambda[g, year]) <- psi_mu +
				inprod(occ_covs[g,1:nlatent ], beta_occ[1:nlatent]) +
				cell_area
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
		background[year] <- inprod(lambda[1:G, year], thin_prob[1:G, year] ) / npo[year]
		zsum[year] <- sum(z[1:G,year])
	}
	#
	# Analyze the (opp)orntunistic PO data
	#  Since there is varying amounts of data per year
	#  we use nested indexing to collect the appropriate parameters
	#  for each PO datapoint.
	#  This means that we have grouped all of the PO data into one long vector.
	#  of length(all_npo)
	# for(opp in 1:all_npo){
	# 	# lambda for PO, nested indexing
	# 	log(lambda_po[opp]) <- psi_mu + psi_season[opp_year[opp]] +
	# 		inprod(occ_covs[po_pixel[opp],1:nlatent ], beta_occ[1:nlatent]) + 
	# 		cell_area
	# 	# thin for PO, nested indexing
	# 	logit(thin_lambda_po[opp]) <- po_mu + po_season[opp_year[opp]] + 
	# 		inprod(po_det_covs[po_pixel[opp], 1:nobs_po], beta_po_det[1:nobs_po])
	# 	# one's trick
	# 	ones[opp] ~ dbern(
	# 		exp(
	# 			log(lambda_po[opp] * thin_lambda_po[opp]) -
	# 				log(background[opp_year[opp]])
	# 		) / CONSTANT
	# 	)
	# }
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
	psi_mu ~ dnorm(-3.2, 1)T( -8.5, 5)
	for(latent in 1:nlatent){
		beta_occ[latent] ~  dnorm(0, 0.75)
	}
	# priors for presence absence observation 
	# temporal random effect
	pa_mu ~ dnorm(0, 0.75)T( -5,5)
	for(obs_pa in 1:nobs_pa){
		beta_pa_det[obs_pa] ~ dnorm(0, 0.75)
	}
	# priors for presence only observation
	# temporal random effect
	po_mu ~ dnorm(0, 0.75)T(-5,5)
	#lambda_po_det ~ dunif(0.001, 20)
	for(obs_po in 1:nobs_po){
		beta_po_det[obs_po]  ~  dnorm(0, 0.75)
	}
}


