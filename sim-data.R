###################################################################################################
## original code by Robert Dorazio
## modified by Vira Koshkina
## Integrated Species Distribution Models: Combining presence-only data and presence-absence data with imperfect detection## Single indepandant covariate for each p and lambda.
## utility functions for the file POPA-functions
## 19/08/2016
###################################################################################################
require(raster)
require(mvtnorm)


source("functions.R")

J.so = 3  #number of repeated surveys
minrecipCondNum = 1e-6
factor = 4
beta.soram = c(log(8000), 0.5)
alpha.soram = c(-1.0,-1.0) #po
alpha.so = c(0,-1.5) #probability of detection for repeated surveys

# Define a rectangular region S
s.xmin = -1
s.xmax =  1
s.ymin = -1
s.ymax =  1

s.area =  (s.xmax-s.xmin)*(s.ymax-s.ymin)

npixels.x = 1000
npixels.y = 1000

s = raster(ncol=npixels.x, nrow=npixels.y, xmn=s.xmin, xmx=s.xmax, ymn=s.ymin, ymx=s.ymax)
s.loc = xyFromCell(s, 1:ncell(s))



# Simulate covariate values over discretization of region S
mu1.x = s.xmin + 0.75*(s.xmax-s.xmin)
mu1.y = s.ymin + 0.40*(s.ymax-s.ymin)
sigma1.x = 0.25*abs(s.xmax-s.xmin)
sigma1.y = 0.50*abs(s.ymax-s.ymin)
rho1.xy = 0.5
mu1 = c(mu1.x, mu1.y)
Sigma1 = matrix(c(sigma1.x^2, rep(rho1.xy*sigma1.x*sigma1.y, 2), sigma1.y^2), ncol=2)

mu2.x = s.xmin + 0.15*(s.xmax-s.xmin)
mu2.y = s.ymin + 0.80*(s.ymax-s.ymin)
sigma2.x = 0.50*abs(s.xmax-s.xmin)
sigma2.y = 0.25*abs(s.ymax-s.ymin)
rho2.xy = -0.4
mu2 = c(mu2.x, mu2.y)
Sigma2 = matrix(c(sigma2.x^2, rep(rho2.xy*sigma2.x*sigma2.y, 2), sigma2.y^2), ncol=2)

xcov = 0.4 * dmvnorm(s.loc, mean=mu1, sigma=Sigma1) + 0.6 * dmvnorm(s.loc, mean=mu2, sigma=Sigma2)
xcov = (xcov - mean(xcov))/sd(xcov)

values(s) = xcov
names(s) = 'x'
s.occupancy=s

mu3.x = s.xmin + 0.25*(s.xmax-s.xmin)
mu3.y = s.ymin + 0.65*(s.ymax-s.ymin)
sigma3.x = 0.25*abs(s.xmax-s.xmin)
sigma3.y = 0.50*abs(s.ymax-s.ymin)
rho3.xy = 0.1
mu3 = c(mu3.x, mu3.y)
Sigma3 = matrix(c(sigma3.x^2, rep(rho3.xy*sigma3.x*sigma3.y, 2), sigma3.y^2), ncol=2)

wcov = dmvnorm(s.loc, mean=mu3, sigma=Sigma3)
wcov = (wcov - mean(wcov))/sd(wcov)


s.detection = raster(s)
values(s.detection) = wcov
## values(temp) = xcov   #   NOTE:  UNCOMMENT THIS LINE TO MAKE COVARIATES x and w IDENTICAL
names(s.detection) = 'w'
s = addLayer(s, s.detection)




# Assign values of model parameters and design matrices


X = cbind(rep(1, ncell(s)), values(s)[,'x'])

W = cbind(rep(1, ncell(s)), values(s)[,'w'])


# ... compute expected limiting density of individuals and true detection probabilities over discretization of region S
temp = raster(s)
values(temp) = exp(X %*% beta.soram)
names(temp) = 'lambda'
s = addLayer(s, temp)

temp = raster(s)
values(temp) = expit(W %*% alpha.soram)
names(temp) = 'pTrue'
s = addLayer(s, temp)

maxlambda = max(values(s)[,'lambda'])

#simulate point-process data --------------------------------------------------------------------

# ... simulate point pattern of individuals over discretization of S
N.hpp = rpois(1, maxlambda*s.area)

ind.hpp = sample(1:ncell(s), size=N.hpp, replace=FALSE)   #  sampling w/o replacement ensures only 1 indiv per pixel
loc.hpp = s.loc[ind.hpp, ]
lambda.hpp = values(s)[,'lambda'][ind.hpp]

ind.ipp = runif(N.hpp, 0,1) <= lambda.hpp/maxlambda
N.ipp = sum(ind.ipp)
loc.ipp = loc.hpp[ind.ipp, ]
plot(s$lambda)
points(loc.ipp)

# ... simulate presence-only data (= detections of individuals as a thinned point process)

pTrue.ipp = values(s)[,'pTrue'][ind.hpp][ind.ipp]
y.ipp = rbinom(length(pTrue.ipp), size=1, prob=pTrue.ipp)


#censoring of the po points --------------------------------------------------------------------------------
#censoring 70 % of the po points
y.ipp = ifelse(runif(length(y.ipp),0,1)<0.7,0,y.ipp)


# ----------------------------------------------------------------------------------------------------------

ind.pb = (y.ipp == 1)
loc.pb = loc.ipp[ind.pb, ]


# ... analyze simulated presence-only data

# ... first establish discrete grid for integrating over region S;
# ... then compute area, average value of covariates, and number of individuals in each element of grid

sgrid = aggregate(s, fact=2, fun=mean)

sgrid.loc = xyFromCell(sgrid, 1:ncell(sgrid))
sgrid.xcov = values(sgrid)[,'x']
sgrid.wcov = values(sgrid)[,'w']


area.back = rep(res(sgrid)[1]*res(sgrid)[2], ncell(sgrid))
X.back = cbind(rep(1, ncell(sgrid)), sgrid.xcov)
W.back = cbind(rep(1, ncell(sgrid)), sgrid.wcov)


X.pb = X[ind.hpp[ind.ipp][ind.pb], ]  # thinned observations
W.pb = W[ind.hpp[ind.ipp][ind.pb], ]  # thinned observations

pb.occupancy=X.pb[,-1]
pb.detection=W.pb[,-1]

#simulate point-count data --------------------------------------------------------------------

# .... first compute abundances within a sample frame
spop = dropLayer(s, c('lambda','pTrue'))
temp = raster(spop)
z = rep(0, ncell(spop))
z[ ind.hpp[ind.ipp] ] = 1
values(temp) = z
names(temp) = 'presence'
spop = addLayer(spop, temp)





# .... form sample frame (= grid of sample units) by aggregating every no. cols and every no. rows in spop
gridfact = c(factor,factor)
sgrid = aggregate(spop, fact=gridfact, fun=mean)

abund = aggregate(subset(spop, 'presence'), fact=gridfact, fun=sum)



# abund=crop(abund, c(xleft,xright,ybottom,ytop))

names(abund) = 'N'
sgrid = addLayer(sgrid, abund)
rm(abund)


#Plot a restricted area on the raster plot
#   par(mfrow=c(1,2))
#   plot(s,1)
#   rect(xleft, ybottom ,xright, ytop)
#   plot(s,2)
#   rect(xleft, ybottom,xright, ytop)

sgrid.loc = xyFromCell(sgrid, 1:ncell(sgrid))
sgrid.xcov = values(sgrid)[,'x']
sgrid.wcov = values(sgrid)[,'w']
sgrid.N = values(sgrid)[,'N']

n.so =200
ind.so = sample(1:ncell(sgrid), size=n.so, replace=FALSE)  # simple random sample w/o replacement
N.so = values(sgrid)[,'N'][ind.so]
area.so = rep(res(sgrid)[1]*res(sgrid)[2], length(N.so))
xcov.so = values(sgrid)[,'x'][ind.so]
wcov.so = values(sgrid)[,'w'][ind.so]
so.occupancy=xcov.so


# .... compute observed counts


W.so = array(dim=c(n.so, J.so, 2))
W.so[,,1] = 1
W.so[,,2] = matrix(rep(wcov.so,J.so), ncol=J.so)
Wmat.so = cbind(rep(1,n.so), wcov.so)
so.detection=W.so[,,2]

pTrue.so = matrix(nrow=n.so, ncol=J.so)
y.so = matrix(nrow=n.so, ncol=J.so)
ocupancy.so=ifelse(N.so>1,1,N.so)
for (j in 1:J.so) {
	pTrue.so[, j] = expit(as.matrix(W.so[,j,], nrow=n.so) %*% alpha.so)
	y.so[, j] = rbinom(n.so, size=ocupancy.so, prob=pTrue.so[, j])
}

#change to PA ----------------------
y.so = ifelse(y.so>1,1,y.so) #simply replace all the values that're larger than one with ones
y.so.pres = y.so[rowSums(y.so)>=1,] #detection/non detection matrix for sites with detection in at least one of the surveys

save(s.occupancy, s.detection, pb.occupancy,pb.detection, y.so, so.occupancy, so.detection, file="data.rda")

