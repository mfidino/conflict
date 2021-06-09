model{
	
	# I'm writing this model in a way that is slightly less comprehensible
	#  But it runs about 20% faster.
	for(g in 1:G){
		for(year in 1:nyear){
		# log-linear predictor intensity across landscape
		#  Need to calculate for each cell because we need
		#  the background estimate for presence only data
		log(lambda[g, year]) <- inprod(X[g,], b[,year]) +
			inprod(occ_covs[g,1:nlatent ], beta_occ[1:nlatent]) +
			cell_area
			
		# logit-linear predictor for thinned Poisson process. 
		#  Similar to lambda, we need to calculate this for each
		#  cell for the background estimate.
		logit(thin_prob[g, year]) <- po_mu + po_season[year] + 
			inprod(po_det_covs[g,1:nobs_po ] , beta_po_det[1:nobs_po])
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
		logit(pa_det_prob[site, year]) <- pa_mu + pa_season[year] +
			inprod(pa_det_covs[pa_pixel[site],1:nobs_pa],beta_pa_det[1:nobs_pa])
		y_pa[site, year] ~ dbin(
			z[pa_pixel[site], year] * pa_det_prob[site, year],
			J[site,year]
		)
		}
	}
	# priors for latent state
	# temporal random effect
	#psi_mu ~ dnorm(-3.2, 0.75)
  # slope terms
	for(latent in 1:nlatent){
			beta_occ[latent] ~ dnorm(0, 0.75)
	}
	# Smoothing spline for occupancy
	b[1,1] ~ dnorm(-3.2,0.75)
	## prior for s(E,N)... 
	K1 <- S1[1:(nspline-1),1:(nspline-1)] * lambda_gam 
	b[2:nspline,1] ~ dmnorm(zero[2:nspline],K1) 
	## smoothing parameter priors CHECK...
	lambda_gam ~ dgamma(.05,.005)
	rho <- log(lambda_gam)
	# yearly variation for spline for spline as random walk
	gam_tau ~ dgamma(1,1)
	gam_sd <- 1 / sqrt(gam_tau)
	for(jj in 1:nspline){
	  for(yr in 2:nyear){
		  b[jj,yr] ~ dnorm(b[jj,yr-1], gam_tau)
	  }
	}
	# priors for presence absence observation 
	# temporal random effect
	pa_mu ~ dnorm(0, 0.75)
	pa_tau_season ~ dgamma(1,1)
	pa_sd_season <- 1 / sqrt(pa_tau_season)
	for(year in 1:nyear){
		pa_season[year] ~ dnorm(0, pa_tau_season)
	}
	for(obs_pa in 1:nobs_pa){
		beta_pa_det[obs_pa] ~ dnorm(0, 0.75)
	}
	# priors for presence only observation
	# temporal random effect
	po_mu ~ dnorm(0, 0.75)
	po_tau_season ~ dgamma(1,1)
	po_sd_season <- 1 / sqrt(po_tau_season)
	#lambda_po_det ~ dunif(0.001, 20)
	for(year in 1:nyear){
		po_season[year] ~ dnorm(0, po_tau_season)
	}
	for(obs_po in 1:nobs_po){
			beta_po_det[obs_po] ~ dnorm(0,0.75)
		}
	}



