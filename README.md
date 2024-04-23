# HEATsims
This repository includes code to reproduce simulation studies from ["Heterogeneity-aware integrative regression for ancestry-specific association studies"](https://arxiv.org/abs/2306.05571). For shorthand, we refer to our estimator as HEAT. 

## Simulation studies 
The simulation studies are designed to be performed on a high-performance computing cluster using slurm. To run the simulations, one needs to create a directory with all the files and directory included in this repository, then execute ```Run_Simulations.sh```. This will initiate 2000 jobs: each performs one replicate of the simulation study by running the ```Protein9_Main.R``` script with a distinct set of model settings and seeds. 

To exactly recreate our results, one would need access to the WHI SNP genotype data (e.g., through dbGaP). For users with access to these data, please email the first author who can provide specific instructions for creating the SNP matrix used in the simulation studies. For those wanting to perform related simulation studies, one need only comment out line 47--89 of "Protein9_Main.R" and create two SNP matrices ``SNP.AA`` and ``SNP.EA``. For example, to create normally distributed predictors with AR(1) covariance structures, use
```
nAA <- 450; nEA <- 900; p <- 500
Sigma.AA <- matrix(0, nrow=p, ncol=p)
Sigma.EA <- matrix(0, nrow=p, ncol=p)
for (kk in 1:p) {
  for (jj in 1:p) {
    Sigma.AA[kk,jj] <- 0.5^abs(jj-kk); Sigma.EA[kk,jj] <- 0.95^abs(jj-kk)
  }
}
eo.AA <- eigen(Sigma.AA); eo.EA <- eigen(Sigma.EA)
SNP.AA <- matrix(rnorm(nAA*p), nrow = nAA) %*% eo.AA$vec %*% diag(eo.AA$val^0.5) %*% t(eo.AA$vec)
SNP.EA <- matrix(rnorm(nEA*p), nrow = nEA) %*% eo.EA$vec %*% diag(eo.EA$val^0.5) %*% t(eo.EA$vec)
```

## Implementing HEAT 
For users simply wanting to implement HEAT, we recommend following the example provided in the ```Example.rmd``` file. At present, the software works with $J = 2$ populations. Please check back for updates after publication. 
