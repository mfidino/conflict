
model{
log(y) ~ dgamma(1,1)
#k <- log(y)
j <- 1 - exp(-exp(y))
}		

