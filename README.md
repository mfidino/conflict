# A repository for:

Fidino, M, Lehrer, E. W., Kay, C. A. M., Yarmey, N., Murray, M. H., Fake, K., Adams, H. C., & Magle, S. B. Combining nuisance wildlife reports with wildlife monitoring data to estimate the probability of human-wildlife conflict relative to a species’ underlying distribution.



![Rought draft of figure 1 from the manuscript. The left plot shows a species distribution, which is highest on the lower left of the plot. The center plot shows a species conflict potential, which is highest on the right side of the plot. The right plot shows where conflict is most likely to occur on the landscape, and is a product of the left and center plots. As such, there is a smaller "hot spot" towards the center of this plot.](./figures/rough_1.png)

This repository stores all of the data and code used to fit the integrated model to the Chicago, Illinois nuisance wildlife complaint data and the camera trap data we collected between 2011 and 2013.

Right now the associated scripts on this repo will simulate and analyze presence-only data with detection / non-detection data from planned surveys. Models to estimate the set parameter values from the simulated data are written in `JAGS`. 

After cloning the repository all you will need to open up is the `simulate_poisson_process.R` script and run through it. The integrated model does take some time to run.
