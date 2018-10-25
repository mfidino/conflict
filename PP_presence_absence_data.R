model{
	
	for(g in 1:G){
		u[g] ~ dpois(lambda[g] * cell_area[g])
		lambda[g] <- exp(inprod(occ_covs[g,], beta_occ))
		psi[g] <- 1 - exp(-lambda[g]*cell_area[g])
		z[g] ~ dbern(psi[g])
		
	}
	
	for(i in 1:nsite){
		logit(p[i]) <- inprod(det_covs[pixel[i],],
													beta_det)
		mu[i] <- z[pixel[i]] * p[i]
		y[i] ~ dbin(mu[i], 4)
	}
	
	beta_occ[1] ~ dlogis(0, 1)
	beta_occ[2] ~ dlogis(0, 1)
	beta_det[1] ~ dlogis(0, 1)
	
}


