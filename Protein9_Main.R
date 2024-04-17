# -------------------------------------------------------------
# Run simple for-loop with 100 random training/testing splits
# -------------------------------------------------------------
namedList <- function(...) {
  L <- list(...)
  snm <- sapply(substitute(list(...)),deparse)[-1]
  if (is.null(nm <- names(L))) nm <- snm
  if (any(nonames <- nm=="")) nm[nonames] <- snm[nonames]
  setNames(L,nm)
}
library(glmnetUtils)
path <- "~/yourpath/yourdirectory/"
library(natural, lib.loc = paste(path, "packages/", sep=""))
source(paste(path, "Functions/penMLE_cvL.R", sep=""))
source(paste(path, "Functions/penMLE_fixedRhos_cvL.R", sep=""))

TPR <- function(x, y){
  sum(x != 0 & y != 0)/sum(x!=0)
}

MS <- function(a){
  sum(a!=0)
}

ModelErr <- function(a, b, Sigma){
  tcrossprod(crossprod(a - b, Sigma), a-b)
}

getResults <- function(beta, estimate, SNP.test, name){
  mse <- sum((beta - estimate)^2)
  tpr <- TPR(beta, estimate)
  ms <- MS(estimate)
  me <- ModelErr(beta, estimate, SNP.test)
  result <- named.list(mse, tpr, ms, me)
  names(result) <- names(paste(result, name, sep=""))
  return(result)
}

# -------------------------------------
# Set protein index
# -------------------------------------
uu <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# ------------------------------------
# Try example here
# ------------------------------------
temp <- readRDS("~/blue/pQTL_Sims/Simulations_R1/Protein9/Protein9.RDS")
SNP.EA <- temp$SNP.EA
SNP.AA <- temp$SNP.AA

# Model 1
keep <- which(apply(SNP.EA, 2, sd)!=0 | apply(SNP.AA, 2, sd)!=0)
SNP.EA <- SNP.EA[,keep]
SNP.AA <- SNP.AA[,keep]
t.AA <- apply(SNP.AA, 2, function(x){min(1 - mean(x)/2, mean(x)/2)})
t.EA <- apply(SNP.EA, 2,  function(x){min(1 - mean(x)/2, mean(x)/2)})
keep <- union(which(t.AA > 0.05), which(t.EA > 0.05))
SNP.EA <- SNP.EA[,keep]
SNP.AA <- SNP.AA[,keep]


SNP_Prune <- function(X, cor.cutoff){
  # sort by MAF
  keep <- NULL
  Xsorted <- X[,sort(pmin(apply(X, 2, mean), 1 - apply(X, 2, mean)), decreasing=TRUE, index=TRUE)$ix]
  temp <- abs(cor(Xsorted[,1], Xsorted[,-1]))
  k <- dim(Xsorted)[2]
  keepnames <- NULL
  while(k > 1){
    keep <- cbind(keep, Xsorted[,1])
    keepnames <- c(keepnames, colnames(Xsorted)[1])
    rm <- which(temp >= cor.cutoff)
    Xsorted <- Xsorted[,-c(1, rm+1), drop=FALSE]
    if(dim(Xsorted)[2] == 0){
      k = 1
    } else {
      temp <- abs(cor(Xsorted[,1], Xsorted[,-1]))
      k <- dim(Xsorted)[2]
      #cat(k, "\n")
    }
  }
  return(list("X" = keep, "snpnames" = keepnames))
}

ea.keep <- SNP_Prune(SNP.EA, 0.95)
aa.keep <- SNP_Prune(SNP.AA, 0.95)
keepnames <- unique(c(ea.keep$snpnames), c(aa.keep$snpnames))
SNP.EA <- SNP.EA[,keepnames]
SNP.AA <- SNP.AA[,keepnames]
n1 <- dim(SNP.AA)[1]
n2 <- dim(SNP.EA)[1]

nreps <- 100
t0 <- expand.grid( 
  "Prop_pQTLsShared" = c(0.6, 0.8, 0.9, 1.0), 
  "sigmaDiff" = c(1.0, 1.25, 1.5, 2, 3),
  "R2" = rep(0.5, each = nreps))

R2 <- t0[uu,3]
Prop_pQTLsShared <- t0[uu,1]
sigmaDiff <- t0[uu,2]
s <- 40
p <- length(keepnames)
aa.sd <- apply(SNP.AA, 2, sd)
ea.sd <- apply(SNP.EA, 2, sd)
beta.AA <- rep(0, p)

set.seed(uu)
preds <- sample(1:p, s)
beta.AA[preds] <- aa.sd[preds]

beta.EA <- rep(0, p)
preds2 <- sample(preds, s*Prop_pQTLsShared)
beta.EA[preds2] <- ea.sd[preds2]
if(Prop_pQTLsShared < 1){
  preds3 <- sample(c(1:p)[-preds2], s - s*Prop_pQTLsShared)
  beta.EA[preds3] <- ea.sd[preds3]
}
beta.EA[which(ea.sd == 0)] <- 0

sign.vec.AA <- sample(c(1, -1), p, replace=TRUE, prob = c(0.8, 0.2))
sign.vec.EA <- sample(c(1, -1), p, replace=TRUE, prob = c(0.8, 0.2))
beta.AA <- beta.AA*sign.vec.AA
beta.EA <- beta.EA*sign.vec.EA 

sigma <- sd(SNP.EA%*%beta.EA)
prot.AA <- -mean(SNP.AA%*%beta.AA) + SNP.AA%*%beta.AA + rnorm(n1, sd = sigma*sigmaDiff)
prot.EA <- -mean(SNP.EA%*%beta.EA) + SNP.EA%*%beta.EA + rnorm(n2, sd = sigma)
SNP.AA.test <- cov(SNP.AA)
SNP.EA.test <- cov(SNP.EA)

me.null <- tcrossprod(crossprod(beta.AA, SNP.AA.test), beta.AA)

# -------------------------------------------
# Initialize results
# -------------------------------------------
Results <- NULL 

# ------------------------------------------
# Combining data (naive approach)
# ------------------------------------------
nfolds <- 10
fold.id <- c(sample(rep(1:nfolds, length=length(prot.EA))), sample(rep(nfolds:1, length=length(prot.AA))))

# ------------------------------------------
# Separate models (SEN) 
# ------------------------------------------
fit.full <- cva.glmnet(y = prot.AA, 
                      x = SNP.AA, nfolds = 10, 
                      alpha = seq(1, 0, len = 11)^3,
                      foldid = fold.id[-(1:length(prot.EA))])
errs <- matrix(Inf, nrow = 100, ncol = 11)
for(kk in 1:11){
  errs[1:length(fit.full$modlist[[kk]]$cvm),kk] <- fit.full$modlist[[kk]]$cvm
}
alpha.ind <- which(errs == min(errs), arr.ind=TRUE)[1,2]
fit.full <- fit.full$modlist[[alpha.ind]]
beta.AA.sep.min <- coef(fit.full, s=fit.full$lambda.min)[-1]
Results <- append(Results, getResults(beta.AA, beta.AA.sep.min, SNP.AA.test, ".AA.sep.min"))

fit.full <- cva.glmnet(y = prot.EA, 
                       x = SNP.EA, nfolds = 10, 
                       alpha = seq(1, 0, len = 11)^3,
                       foldid = fold.id[(1:length(prot.EA))])
errs <- matrix(Inf, nrow = 100, ncol = 11)
for(kk in 1:11){
  errs[1:length(fit.full$modlist[[kk]]$cvm),kk] <- fit.full$modlist[[kk]]$cvm
}
alpha.ind <- which(errs == min(errs), arr.ind=TRUE)[1,2]
fit.full <- fit.full$modlist[[alpha.ind]]
beta.EA.sep.min <- coef(fit.full, s=fit.full$lambda.min)[-1]
Results <- append(Results, getResults(beta.EA, beta.EA.sep.min, SNP.EA.test, ".EA.sep.min"))


# ----------------------------------------
# Combining data (naive approach)
# ----------------------------------------
fit.full <- cva.glmnet(y = c(prot.EA, prot.AA), 
                      x = rbind(cbind(1, SNP.EA), cbind(0, SNP.AA)),
                      alpha = seq(1, 0, len = 11)^3,
                      penalty.factor = c(0, rep(1, dim(SNP.EA)[2])),
                      nfolds = 10, foldid = fold.id)

errs <- matrix(Inf, nrow = 100, ncol = 11)
for(kk in 1:11){
  errs[1:length(fit.full$modlist[[kk]]$cvm),kk] <- fit.full$modlist[[kk]]$cvm
}
alpha.ind <- which(errs == min(errs), arr.ind=TRUE)[1,2]
fit.full <- fit.full$modlist[[alpha.ind]]
beta.AA.both.min <- coef(fit.full, s=fit.full$lambda.min)[-c(1,2)]
Results <- append(Results, getResults(beta.AA, beta.AA.both.min, SNP.AA.test, ".AA.both.min"))

# --------------------------------------------------
# Our estimator 
# -------------------------------------------------
ptm <- proc.time()
nfolds <- 10
fold.id.AA <- sample(rep(nfolds:1, length=length(prot.AA)))
fold.id.EA <- sample(rep(nfolds:1, length=length(prot.EA)))
temp <- penMLE_CV(prot.AA, prot.EA, SNP.AA, SNP.EA, 
                  nlambda = 20, delta = 0.01, ngamma = 20,
                  nfolds = 10, fold.id.AA = fold.id.AA, 
                  fold.id.EA = fold.id.EA)
getTime <- proc.time() - ptm
times.ours.min <- getTime[3]

# using sum of squared residuals as tuning criterion 
errs <- apply(temp$errs.AA, c(1,2), sum) + apply(temp$errs.EA, c(1,2), sum) 
t0 <- which(errs == min(errs), arr.ind=TRUE)
getInd <- t0[1,2]
full.fit <- temp$full.fit[[t0[1,1]]]
full.sd <- apply(rbind(SNP.AA, SNP.EA), 2, sd)
beta.AA.ours.min <- (full.fit$phiA[,getInd]/full.sd)/full.fit$rhoA[getInd]
Results <- append(Results, getResults(beta.AA, beta.AA.ours.min, SNP.AA.test, ".AA.ours.min"))
beta.EA.ours.min <- (full.fit$phiE[,getInd]/full.sd)/full.fit$rhoE[getInd]
Results <- append(Results, getResults(beta.EA, beta.EA.ours.min, SNP.EA.test, ".EA.ours.min"))

# using validation likelihood as tuning criterion 
L.errs <- apply(temp$L.errs.AA, c(1,2), sum) + apply(temp$L.errs.EA, c(1,2), sum) 
t0 <- which(L.errs == min(L.errs), arr.ind=TRUE)
getInd <- t0[1,2]
full.fit <- temp$full.fit[[t0[1,1]]]
full.sd <- apply(rbind(SNP.AA, SNP.EA), 2, sd)
beta.AA.ours.Lmin <- (full.fit$phiA[,getInd]/full.sd)/full.fit$rhoA[getInd]
Results <- append(Results, getResults(beta.AA, beta.AA.ours.Lmin, SNP.AA.test, ".AA.ours.Lmin"))

beta.EA.ours.Lmin <- (full.fit$phiE[,getInd]/full.sd)/full.fit$rhoE[getInd]
Results <- append(Results, getResults(beta.EA, beta.EA.ours.Lmin, SNP.EA.test, ".EA.ours.Lmin"))


# ---------------------------------------
# reHEAT 
# ---------------------------------------
temp2 <- penMLE_fixedrhos_CV(prot.AA, prot.EA, SNP.AA, SNP.EA, 
                            nlambda = 20, delta = 0.01, ngamma = 20,
                            nfolds = 10, fold.id.AA = fold.id.AA, 
                            fold.id.EA = fold.id.EA, 
                            rhoA = 1, rhoE = 1)

errs2 <- apply(temp2$errs.AA, c(1,2), sum) + apply(temp2$errs.EA, c(1,2), sum)
t0 <- which(errs2 == min(errs2), arr.ind=TRUE)
getInd <- t0[1,2]
full.fit <- temp2$full.fit[[t0[1,1]]]
full.sd <- apply(rbind(SNP.AA, SNP.EA), 2, sd)
beta.AA.fixed1.min <- (full.fit$phiA[,getInd]/full.sd)
Results <- append(Results, getResults(beta.AA, beta.AA.fixed1.min, SNP.AA.test, ".AA.fixed1.min"))

beta.EA.fixed1.min <- (full.fit$phiE[,getInd]/full.sd)
Results <- append(Results, getResults(beta.EA, beta.EA.fixed1.min, SNP.EA.test, ".EA.fixed1.min"))



# ---------------------------------------
# Oracle-heat
# ---------------------------------------
temp3 <- penMLE_fixedrhos_CV(prot.AA, prot.EA, SNP.AA, SNP.EA, 
                            nlambda = 20, delta = 0.01, ngamma = 20,
                            nfolds = 10, fold.id.AA = fold.id.AA, 
                            fold.id.EA = fold.id.EA, 
                            rhoA = 1/(sigma*sigmaDiff), 
                            rhoE = 1/sigma)

errs3 <- apply(temp3$errs.AA, c(1,2), sum) + apply(temp3$errs.EA, c(1,2), sum)
t0 <- which(errs3 == min(errs3), arr.ind=TRUE)
getInd <- t0[1,2]
full.fit <- temp3$full.fit[[t0[1,1]]]
full.sd <- apply(rbind(SNP.AA, SNP.EA), 2, sd)
beta.AA.fixedOracle.min <- (full.fit$phiA[,getInd]/full.sd)/full.fit$rhoA
Results <- append(Results, getResults(beta.AA, beta.AA.fixedOracle.min, SNP.AA.test, ".AA.fixedOracle.min"))

beta.EA.fixedOracle.min <- (full.fit$phiE[,getInd]/full.sd)/full.fit$rhoE
Results <- append(Results, getResults(beta.EA, beta.EA.fixedOracle.min, SNP.EA.test, ".EA.fixedOracle.min"))


errs3 <- apply(temp3$L.errs.AA, c(1,2), sum) + apply(temp3$L.errs.EA, c(1,2), sum)
t0 <- which(errs3 == min(errs3), arr.ind=TRUE)
getInd <- t0[1,2]
full.fit <- temp3$full.fit[[t0[1,1]]]
full.sd <- apply(rbind(SNP.AA, SNP.EA), 2, sd)
beta.AA.fixedOracle.Lmin <- (full.fit$phiA[,getInd]/full.sd)/full.fit$rhoA
Results <- append(Results, getResults(beta.AA, beta.AA.fixedOracle.Lmin, SNP.AA.test, ".AA.fixedOracle.Lmin"))

beta.EA.fixedOracle.Lmin <- (full.fit$phiE[,getInd]/full.sd)/full.fit$rhoE
Results <- append(Results, getResults(beta.EA, beta.EA.fixedOracle.Lmin, SNP.EA.test, ".EA.fixedOracle.Lmin"))

# -------------------------------------------------------------------
# Natural lasso estimators of variances, then our method (HEAT-App)
# -------------------------------------------------------------------
ptm <- proc.time()
sigma.AA <- nlasso_cv(x = SNP.AA, y = prot.AA, nfold=10)$sig_obj
sigma.EA <- nlasso_cv(x = SNP.EA[,which(apply(SNP.EA, 2, sd)!=0)], y = prot.EA, nfold=10)$sig_obj
temp4 <- penMLE_fixedrhos_CV(prot.AA, prot.EA, SNP.AA, SNP.EA, 
                             nlambda = 20, delta = 0.01, ngamma = 20,
                             nfolds = 10, fold.id.AA = fold.id.AA, 
                             fold.id.EA = fold.id.EA, 
                             rhoA = 1/sigma.AA, 
                             rhoE = 1/sigma.EA)

getTime <- proc.time() - ptm
times.oneStep.min <- getTime[3]

errs4 <- apply(temp4$errs.AA, c(1,2), sum) + apply(temp4$errs.EA, c(1,2), sum)
t0 <- which(errs4 == min(errs4), arr.ind=TRUE)
getInd <- t0[1,2]
full.fit <- temp4$full.fit[[t0[1,1]]]
full.sd <- apply(rbind(SNP.AA, SNP.EA), 2, sd)

beta.AA.oneStep.min <- (full.fit$phiA[,getInd]/full.sd)/full.fit$rhoA
Results <- append(Results, getResults(beta.AA, beta.AA.oneStep.min, SNP.AA.test, ".AA.oneStep.min"))


beta.EA.oneStep.min <- (full.fit$phiE[,getInd]/full.sd)/full.fit$rhoE
Results <- append(Results, getResults(beta.EA, beta.EA.oneStep.min, SNP.EA.test, ".EA.oneStep.min"))


errs4 <- apply(temp4$L.errs.AA, c(1,2), sum) + apply(temp4$L.errs.EA, c(1,2), sum)
t0 <- which(errs4 == min(errs4), arr.ind=TRUE)
getInd <- t0[1,2]
full.fit <- temp4$full.fit[[t0[1,1]]]
full.sd <- apply(rbind(SNP.AA, SNP.EA), 2, sd)

beta.AA.oneStep.Lmin <- (full.fit$phiA[,getInd]/full.sd)/full.fit$rhoA
Results <- append(Results, getResults(beta.AA, beta.AA.oneStep.min, SNP.AA.test, ".AA.oneStep.Lmin"))


beta.EA.oneStep.Lmin <- (full.fit$phiE[,getInd]/full.sd)/full.fit$rhoE
Results <- append(Results, getResults(beta.EA, beta.EA.oneStep.min, SNP.EA.test, ".EA.oneStep.Lmin"))


# --------------------------------------------
# Save results + model information
# --------------------------------------------
Results <- append(Results, namedList(
  beta.AA, 
  beta.EA, 
  sigma, 
  sigmaDiff, 
  me.null,
  sigma.AA,
  sigma.EA, 
  times.oneStep.min,
  times.ours.min))


# --------------------------------------------
# Save results + model information
# --------------------------------------------
saveRDS(Results, file = paste(path, "Protein9/Protein9_Result/Result", uu,".RDS", sep=""))

q("no")
