model{
	
	for(g in 1:G){
		u[g] ~ dpois(lambda[g])
		lambda[g] <- exp(inprod(occ_covs[g,], beta_occ)) * cell_area[g]
		psi[g] <- 1 - exp(-lambda[g])
		z[g] ~ dbern(psi[g])
		
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
	
}


