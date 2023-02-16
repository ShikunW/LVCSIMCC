
library(foreach)
library(Matrix)

Cost_design <- function(Y, X, time, surv, id, delta, 
                        ndx0=5, deg0=2, ndx=0, deg=1, censor=F,
                        surv.out, family = 'normal'){
  ## B-splines basis ordinate (x)
  tl <- sl <- 0
  tr <- sr <- 1
  nbeta0 <- ndx0 + deg0 + 1
  nbeta <- ndx + deg + 1
  p <- ncol(X)
  nparam <- nbeta0+p
  idx.sts <- which(delta == 1)
  idx.lts <- which(delta == 0 & surv == 1)
  
  Bt0 <- Cost_pbase(time / surv, tl, tr, ndx0, deg0)
  Bs0 <- Cost_pbase(surv, sl, sr, ndx0, deg0)
  B0.sts <- kronecker(matrix(1, ncol = nbeta0), Bt0[idx.sts,]) *
    kronecker(Bs0[idx.sts,], matrix(1, ncol = nbeta0))
  B0.lts <- Bt0[idx.lts,]
  
  Bt <- Cost_pbase(time / surv, tl, tr, ndx, deg)
  Bs <- Cost_pbase(surv, sl, sr, ndx, deg)
  B.sts <- kronecker(matrix(1, ncol = nbeta), Bt[idx.sts,]) *
    kronecker(Bs[idx.sts,], matrix(1, ncol = nbeta))
  B.lts <- Bt[idx.lts,]
  
  # combine results for sts and lts
  s.unc <- surv[c(idx.sts, idx.lts)]
  t.unc <- time[c(idx.sts, idx.lts)]
  delta.unc <- delta[c(idx.sts, idx.lts)]
  id.unc <- id[c(idx.sts, idx.lts)]
  X.unc <- X[c(idx.sts, idx.lts),]
  Y.unc <- Y[c(idx.sts, idx.lts)]
  B0.unc <- as.matrix(bdiag(B0.sts, B0.lts))
  B.unc <- as.matrix(bdiag(B.sts, B.lts))
  
  # calculate conditional survival probability weights
  # collect data for censored patients
  if(!censor){
    time <- t.unc
    surv <- s.unc
    Y <- Y.unc
    X <- X.unc
    id <- id.unc
    delta <- delta.unc
    B <- list(B.unc,c(0,0))
    B0 <- list(B0.unc,c(0,0))
    dB <- 0; dB0 <- 0
  }else{
    # patients censored prior to tau
    idx.cen <- which(delta == 0 & surv < 1)
    t.cen <- time[idx.cen]
    s.cen <- surv[idx.cen]
    Y.cen <- Y[idx.cen]
    X.cen <- X[idx.cen,]
    id.cen <- id[idx.cen]
    delta.cen <- delta[idx.cen]
    
    CondSurv <- surv.out$surv#surv.out$cond_surv$surv0
    CondProb <- surv.out$prob#surv.out$cond_surv$prob
    dprob.cond <- surv.out$cond_surv$dprob
    # uu = match(1:n,id.unique)
    # id.cen1 <- match(id.cen,uu)
    ncen <- length(idx.cen)
    # get design matrix for each subject i=1,...,n; for survival time m=1,...,11
    
    # the following code may change if link function is not identity
    nsurv <- ncol(CondSurv)
    B0.cen.sts1 <- foreach(j = 1:(nsurv-1)) %do% {
      Bt0.cen <- Cost_pbase(t.cen / CondSurv[id.cen, j], tl, tr, ndx0, deg0)
      Bs0.cen <- Cost_pbase(CondSurv[id.cen, j], sl, sr, ndx0, deg0)
      kronecker(matrix(1, ncol = nbeta0), Bt0.cen) *
        kronecker(Bs0.cen, matrix(1, ncol = nbeta0))
    }
    Bt0.cen <- Cost_pbase(t.cen, tl, tr, ndx0, deg0)
    B0.cen.lts1 <- Bt0.cen
    B.cen.sts1 <- foreach(j = 1:(nsurv-1)) %do% {
      Bt.cen <- Cost_pbase(t.cen / CondSurv[id.cen, j], tl, tr, ndx, deg)
      Bs.cen <- Cost_pbase(CondSurv[id.cen, j], sl, sr, ndx, deg)
      kronecker(matrix(1, ncol = nbeta), Bt.cen) *
        kronecker(Bs.cen, matrix(1, ncol = nbeta))
    }
    Bt.cen <- Cost_pbase(t.cen, tl, tr, ndx, deg)
    B.cen.lts1 <- Bt.cen
    
    if(family == 'normal'){
      B0.cen.sts <- foreach(j = 1:(nsurv-1), .combine = '+') %do% {
        CondProb[id.cen, j] * B0.cen.sts1[[j]]
      }
      B0.cen.lts <- CondProb[id.cen, nsurv] * B0.cen.lts1
      B0.cen <- cbind(B0.cen.sts, B0.cen.lts)
      B0 <- rbind(B0.unc, B0.cen)
      B.cen.sts <- foreach(j = 1:(nsurv-1), .combine = '+') %do% {
        CondProb[id.cen, j] * B.cen.sts1[[j]]
      }
      B.cen.lts <- CondProb[id.cen, nsurv] * B.cen.lts1
      B.cen <- cbind(B.cen.sts, B.cen.lts)
      
      B <- rbind(B.unc, B.cen)
    }else if(family %in% c('poisson','binomial')){
      B0.cen <- list(B0.cen.sts, B0.cen.lts)
      B0 <- list(B0.unc, B0.cen)
      B.cen <- list(B.cen.sts, B.cen.lts)
      B <- list(B.unc, B.cen)
    }
    
    time <- c(t.unc, t.cen)
    surv <- c(s.unc, s.cen)
    Y <- c(Y.unc, Y.cen)
    X <- rbind(X.unc, X.cen)
    id <- c(id.unc, id.cen)
    delta <- c(delta.unc, delta.cen)
  }
  out <- list(Y=Y, X=X, time=time, surv=surv,
              id=id, delta=delta, 
              B0=B0, B=B)
  return(out)
}
