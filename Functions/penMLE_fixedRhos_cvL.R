# ------------------------------------------------------------------------------------
# Main function: fits penalized MLE with rhos fixed for fixed tuning parameter pair
# ------------------------------------------------------------------------------------
penMLE_fixedrhos <- function(prot.AA, prot.EA, X.AA, X.EA, phiA, 
                           phiE, rhoA, rhoE, lambda, gamma, alpha0, w){
  
  gradFunction <- function(prot.AA, prot.EA, X.AA, X.EA, X.phiA, X.phiE, rhoA, rhoE){
    n <- length(prot.AA) + length(prot.EA)
    2*cbind(crossprod(X.AA, 
                      X.phiA - prot.AA*rhoA)/n, 
            crossprod(X.EA, 
                      X.phiE - prot.EA*rhoE)/n)
  }
  
  fFunction <- function(prot.AA, prot.EA, X.phiA, X.phiE, rhoA, rhoE){
    n <- length(prot.AA) + length(prot.EA)
    out <- sum((rhoA*prot.AA - X.phiA)^2)/n + 
      sum((rhoE*prot.EA - X.phiE)^2)/n
    return(out)
  }
  
  proxCompute <- function(input, tau, gamma, w){
    t0 <- matrix(0, nrow=nrow(input), ncol=ncol(input))
    input <- pmax(abs(input) - gamma, 0)*sign(input)
    eucRows <- sqrt(rowSums(input^2))
    t1 <- which(eucRows >= w*tau)
    if(length(t1) > 0){
      for(ll in t1){
        t0[ll, ] <- (1 - (tau*w[ll])/sqrt(sum(input[ll,]^2)))*input[ll,]
      } 
    }
    return(t0)
  }
  
  obj.prev <- fFunction(prot.AA, prot.EA, crossprod(t(X.AA), phiA), crossprod(t(X.EA), phiE), rhoA, rhoE) + lambda*sum(sqrt(rowSums(cbind(phiA, phiE)^2)))
  obj.orig <- obj.prev
  alpha0 <- 1
  alpha <- alpha0
  n <- length(prot.AA) + length(prot.EA)
  
  for(k in 1:1e3){
    
    inputs <- which(phiA!=0)
    t1 <- crossprod(t(X.AA[,inputs]),phiA[inputs])
    inputs <- which(phiE!=0)
    t2 <- crossprod(t(X.EA[,inputs]),phiE[inputs])

    tempGrad <-  gradFunction(prot.AA, prot.EA, X.AA, X.EA, t1, t2, rhoA, rhoE)
    linesearch <- TRUE
    tempPrev <- fFunction(prot.AA, prot.EA, t1, t2, rhoA, rhoE)
    
    while(linesearch){
      searchPoint <- cbind(phiA, phiE) - alpha*tempGrad
      phi.temp <- proxCompute(searchPoint, alpha*lambda, alpha*gamma, w)
      
      tempLik <- fFunction(prot.AA, prot.EA, crossprod(t(X.AA), phi.temp[,1]), crossprod(t(X.EA), phi.temp[,2]), rhoA, rhoE)
      tempOld <- (tempPrev + sum(tempGrad*(phi.temp - cbind(phiA, phiE))) + (1/(2*alpha))*sum((phi.temp - cbind(phiA, phiE))^2))
      t1 <- (tempLik <= tempOld)
      if(t1){
        phiAprev <- phiA
        phiEprev <- phiE
        phiA <- phi.temp[,1]
        phiE <- phi.temp[,2]
        linesearch <- FALSE
      } else {
        alpha <- alpha*0.5
      }
    }
    
    tempLik.adj <- tempLik - 2*length(prot.AA)*log(rhoA)/n - 2*length(prot.EA)*log(rhoE)/n
    obj.new <- tempLik.adj + lambda*sum(w*sqrt(rowSums(cbind(phiA, phiE)^2))) + gamma*(sum(abs(phiA)) + sum(abs(phiE)))
    if(k > 5){
      if(abs(obj.new - obj.prev) < 1e-8*abs(obj.prev)){
        break
      }
    }
    obj.prev <- obj.new
  }
  
  return(list("phiA" = phiA,
              "phiE" = phiE, 
              "rhoA" = rhoA, 
              "rhoE" = rhoE))
}

# ----------------------------------------------------------------------------------
# Function for fitting the entire solution path
# stop.index is the number of tuning parameters that must be computed along the path
# if CV errors are decreasing for many consecutive iterations, we terminate
# ------------------------------------------------------------------------------------
penMLE_fixedrhos_path <- function(prot.AA, prot.EA, X.AA, X.EA, lambda, gamma, rhoA, rhoE,
                                prot.AA.test = NULL, prot.EA.test = NULL, 
                                X.AA.test = NULL, X.EA.test = NULL, w, stop.index = 7){
  
  prot.AA.stand <- prot.AA - mean(prot.AA)
  prot.EA.stand <- prot.EA - mean(prot.EA)
  full.sd <- apply(rbind(X.AA, X.EA), 2, sd)
  X.AA.stand <- (X.AA - rep(1, dim(X.AA)[1])%*%t(apply(X.AA, 2, mean)))/(rep(1, dim(X.AA)[1])%*%t(full.sd))
  EA.sd <- apply(X.EA, 2, sd)
  X.EA.stand <- (X.EA - rep(1, dim(X.EA)[1])%*%t(apply(X.EA, 2, mean)))/(rep(1, dim(X.EA)[1])%*%t(full.sd))
  if(any(EA.sd == 0)){
    X.EA.stand[,which(EA.sd == 0)] <- 0
  }                                                                 
  
  phiA <- matrix(0, nrow=dim(X.AA)[2], ncol=length(gamma))
  phiE <- matrix(0, nrow=dim(X.AA)[2], ncol=length(gamma))
  phiA.temp <- phiA[,1]
  phiE.temp <- phiE[,1]
  
  if(is.null(prot.AA.test)){
    for(kk in 1:length(gamma)){
      fit.temp <- penMLE_fixedrhos(prot.AA.stand, prot.EA.stand, X.AA.stand, X.EA.stand, 
                                 phiA = phiA.temp, phiE = phiE.temp, rhoA = rhoA, rhoE = rhoE,
                                 lambda = lambda, gamma = gamma[kk], alpha0 = NULL, w = w)
      phiA[,kk] <- fit.temp$phiA
      phiE[,kk] <- fit.temp$phiE
      phiA.temp <- fit.temp$phiA
      phiE.temp <- fit.temp$phiE
    }
  }

  
  if(!is.null(prot.AA.test)){
    
    full.mean <- apply(rbind(X.AA, X.EA), 2, mean)
    full.sd <- apply(rbind(X.AA, X.EA), 2, sd)
    X.AA.test.stand <- (X.AA.test - rep(1, dim(X.AA.test)[1])%*%t(apply(X.AA, 2, mean)))/(rep(1, dim(X.AA.test)[1])%*%t(full.sd))
    EA.sd <- apply(X.EA, 2, sd)
    X.EA.test.stand <- (X.EA.test - rep(1, dim(X.EA.test)[1])%*%t(apply(X.EA, 2, mean)))/(rep(1, dim(X.EA.test)[1])%*%t(full.sd))
    if(any(EA.sd == 0)){
      X.EA.test.stand[,which(EA.sd == 0)] <- 0
    }
    err.AA <- rep(Inf, length(gamma))
    err.EA <- rep(Inf, length(gamma))
    L.err.AA <- rep(Inf, length(gamma))
    L.err.EA <- rep(Inf, length(gamma))
    for(kk in 1:length(gamma)){
        fit.temp <- penMLE_fixedrhos(prot.AA.stand, prot.EA.stand, X.AA.stand, X.EA.stand, 
                                   phiA = phiA.temp, phiE = phiE.temp, rhoA = rhoA, rhoE = rhoE,
                                   lambda = lambda, gamma = gamma[kk], alpha0=alpha0, w=w)
        phiA[,kk] <- fit.temp$phiA
        phiE[,kk] <- fit.temp$phiE
        phiA.temp <- fit.temp$phiA
        phiE.temp <- fit.temp$phiE
        err.AA[kk] <- sum((prot.AA.test  - mean(prot.AA) - X.AA.test.stand%*%phiA[,kk]/rhoA)^2)
        err.EA[kk] <- sum((prot.EA.test  - mean(prot.EA) - X.EA.test.stand%*%phiE[,kk]/rhoE)^2)
        L.err.AA[kk] <- err.AA[kk]*(rhoA^2) - length(prot.AA.test)*log(rhoA^2)
        L.err.EA[kk] <- err.EA[kk]*(rhoE^2) - length(prot.EA.test)*log(rhoE^2)
        if(kk > stop.index){
          if(all(err.AA[kk:(kk-stop.index-1)] > err.AA[(kk-1):(kk-stop.index)]) & all(err.EA[kk:(kk-stop.index-1)] > err.EA[(kk-1):(kk-stop.index)]) &
            all(L.err.AA[kk:(kk-stop.index-1)] > L.err.AA[(kk-1):(kk-stop.index)]) & all(L.err.EA[kk:(kk-stop.index-1)] > L.err.EA[(kk-1):(kk-stop.index)])){
            break
          }
        }
      }
    return(list("err.AA" = err.AA, 
                "err.EA" = err.EA,
                "L.err.AA" = L.err.AA, 
                "L.err.EA" = L.err.EA))
  } else {
    return(list("phiA" = phiA, "phiE" = phiE, "rhoA" = rhoA, "rhoE"= rhoE))
  }
  
}

# --------------------------------------------
#
# -------------------------------------------

penMLE_fixedrhos_CV <- function(prot.AA, prot.EA, X.AA, X.EA, delta = 0.1, nlambda = 10, 
                              ngamma = 20, nfolds = 5, fold.id.AA, fold.id.EA,
                              rhoA, rhoE){
  
  
  gradFunction <- function(prot.AA, prot.EA, X.AA, X.EA, X.phiA, X.phiE, rhoA, rhoE){
    n <- length(prot.AA) + length(prot.EA)
    2*cbind(crossprod(X.AA, 
                      X.phiA - prot.AA*rhoA)/n, 
            crossprod(X.EA, 
                      X.phiE - prot.EA*rhoE)/n)
  }
  
  prot.AA.stand <- prot.AA - mean(prot.AA)
  prot.EA.stand <- prot.EA - mean(prot.EA)
  full.mean <- apply(rbind(X.AA, X.EA), 2, mean)
  full.sd <- apply(rbind(X.AA, X.EA), 2, sd)
  X.AA.stand <- (X.AA - rep(1, dim(X.AA)[1])%*%t(apply(X.AA, 2, mean)))/(rep(1, dim(X.AA)[1])%*%t(full.sd))
  EA.sd <- apply(X.EA, 2, sd)
  X.EA.stand <- (X.EA - rep(1, dim(X.EA)[1])%*%t(apply(X.EA, 2, mean)))/(rep(1, dim(X.EA)[1])%*%t(full.sd))
  if(any(EA.sd == 0)){
    X.EA.stand[,which(EA.sd == 0)] <- 0
  }
  w <- rep(sqrt(2), dim(X.EA)[2])
  if(any(EA.sd == 0)){
    X.EA.stand[,which(EA.sd == 0)] <- 0
    w[which(EA.sd == 0)] <- 1
  }
  w <- w/max(w)

  tempGrad <- gradFunction(prot.AA.stand, prot.EA.stand, X.AA.stand, X.EA.stand, rep(0, dim(X.AA.stand)[1]), rep(0, dim(X.EA.stand)[1]), rhoA, rhoE)
  t0 <- max(sqrt(rowSums(tempGrad^2)))
  lambda <- 10^seq(log10(0.9*max(abs(t0))), log10(delta*max(abs(tempGrad))), length=nlambda)
  gamma.mat <- matrix(0, nrow=ngamma, ncol=nlambda)
  gamma.try <- 10^seq(-8, 8, length=100)
  for(ll in 1:length(lambda)){
    for(kk in 1:length(gamma.try)){
      t0 <- pmax(abs(tempGrad) - gamma.try[kk], 0)*sign(tempGrad)
      if(max(sqrt(rowSums(t0^2))) < lambda[ll]){
        break
      }
    }
    gamma.mat[,ll] <- 10^seq(log10(gamma.try[kk]), log10(10^(-5)*gamma.try[kk]), length=ngamma)
  }
  
 if(!is.null(nfolds)){

    errs.AA <- array(0, dim=c(nlambda, ngamma, nfolds))
    errs.EA <- array(0, dim=c(nlambda, ngamma, nfolds))
    L.errs.AA <- array(0, dim=c(nlambda, ngamma, nfolds))
    L.errs.EA <- array(0, dim=c(nlambda, ngamma, nfolds))
    t0 <- matrix(0, nrow=nlambda, ncol=ngamma)
    for(kk in 1:nfolds){
      for(ll in 1:nlambda){
        cv.fit <- penMLE_fixedrhos_path(prot.AA = prot.AA[-which(fold.id.AA == kk)], 
                                      prot.EA= prot.EA[-which(fold.id.EA == kk)], 
                                      X.AA = X.AA[-which(fold.id.AA == kk),], 
                                      X.EA = X.EA[-which(fold.id.EA == kk),], 
                                      lambda = lambda[ll], gamma = gamma.mat[which(t0[ll,]!=Inf),ll], 
                                      prot.AA.test = prot.AA[which(fold.id.AA == kk)], 
                                      prot.EA.test = prot.EA[which(fold.id.EA == kk)], 
                                      X.AA.test = X.AA[which(fold.id.AA == kk),], 
                                      X.EA.test = X.EA[which(fold.id.EA == kk),],
                                      rhoA = rhoA,
                                      rhoE = rhoE, w=w)
        errs.AA[ll,,kk] <- cv.fit$err.AA
        errs.EA[ll,,kk] <- cv.fit$err.EA
        L.errs.AA[ll,,kk] <- cv.fit$L.err.AA
        L.errs.EA[ll,,kk] <- cv.fit$L.err.EA
      }
      cat("Through fold ", kk, "\n")
    }

    full.fit <- list()
    for(ll in 1:nlambda){
      full.fit[[ll]] <- penMLE_fixedrhos_path(prot.AA = prot.AA, 
                                      prot.EA = prot.EA, 
                                      X.AA = X.AA, 
                                      X.EA = X.EA, 
                                      lambda = lambda[ll], 
                                      gamma = gamma.mat[which(t0[ll,]!=Inf),ll],
                                      rhoA = rhoA,
                                      rhoE = rhoE, w=w)
    }
    return(list("errs.AA" = errs.AA, "errs.EA" = errs.EA,
    "L.errs.AA" = L.errs.AA, "L.errs.EA" = L.errs.EA, "full.fit" = full.fit,
      "fold.id.AA" = fold.id.AA, "fold.id.EA" = fold.id.EA, 
      "lambda" = lambda, "gamma.mat" = gamma.mat))

  } else {

    full.fit <- list()
    for(ll in 1:nlambda){
      full.fit[[ll]] <- penMLE_fixedrhos_path(prot.AA = prot.AA, 
                                      prot.EA = prot.EA, 
                                      X.AA = X.AA, 
                                      X.EA = X.EA, 
                                      lambda = lambda[ll], 
                                      gamma = gamma.mat[,ll],
                                      rhoA = rhoA,
                                      rhoE = rhoE, w=w)
    }
    return(list("full.fit" = full.fit, 
      "lambda" = lambda, "gamma.mat" = gamma.mat)) 
  }
}

