
model{
y ~ dgamma(1,1)
k <- log(y)
j <- 1 - exp(-exp(y))
}		

 [1] 0.0 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1.0
