model{
	# latent state
	for(g in 1:G){
		for(year in 1:nyear){
		u[g, year] ~ dpois(lambda[g, year])
		lambda[g, year] <- exp(inprod(occ_covs[g,], beta_occ[ ,year])) * cell_area[g]
		psi[g, year] <- 1 - exp(-lambda[g, year])
		z[g, year] ~ dbern(psi[g, year])
		}
	}
	# observation model
	for(site in 1:npa){
		for(year in 1:nyear){
		logit(pa_det_prob[site, year]) <- inprod(pa_det_covs[pa_pixel[site],],
																			 beta_pa_det[ ,year])
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
  # priors for observation 
	#  currently does not vary by year
	for(obs in 1:nobs){
		beta_observation[obs] ~ dlogis(0, 1)
		for(year in 1:nyear){
			beta_pa_det[obs, year] <- beta_observation[obs]
	}
}


