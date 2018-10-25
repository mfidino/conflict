model{
	
	for(g in 1:G){
		probs[g] <- exp(inprod(occ_covs[g,], beta_occ))
	}
	
	for(i in 1:(nind + nz)){
		z[i] ~ dbern(psi)
		pixel[i] ~ dcat(probs[])
	  logit(p[i]) <- inprod(det_covs[pixel[i],],
	  											beta_det)
	  mu[i] <- z[i] * p[i]
	  y[i] ~ dbin(mu[i], 4)
		}
	
	beta_occ[1] ~ dlogis(0, 1)
	beta_occ[2] ~ dlogis(0, 1)
	beta_det[1] ~ dlogis(0, 1)
	psi ~ dbeta(1,1)
	
}


