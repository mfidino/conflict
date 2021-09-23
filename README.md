# A repository for:

Fidino, M, Lehrer, E. W., Kay, C. A. M., Yarmey, N., Murray, M. H., Fake, K., Adams, H. C., & Magle, S. B. Combining nuisance wildlife reports with wildlife monitoring data to estimate the probability of human-wildlife conflict relative to a species’ underlying distribution.

This repository stores all of the data and code used to fit the integrated model to the Chicago, Illinois nuisance wildlife complaint data and the camera trap data we collected between 2011 and 2013.

![longshot](./figures/figure_1.pdf)


Right now the associated scripts on this repo will simulate and analyze presence-only data with detection / non-detection data from planned surveys. Models to estimate the set parameter values from the simulated data are written in `JAGS`. 

After cloning the repository all you will need to open up is the `simulate_poisson_process.R` script and run through it. The integrated model does take some time to run.
