library(mgcv)
library(MASS)
library(Rcpp)
library(RcppArmadillo)
sourceCpp('paraUpdate.cpp')


MAX.IT <- 30
CORE_NUM <- parallel::detectCores()-1

Cost_estimate <- function(Y, X, B0, B, id, delta, surv, CondProb,
                          ndx, deg, ndx0, deg0, ndx1, deg1, 
                          param.init=NULL, lambda, correlation, 
                          varfunc, censor=F, family = 'normal',constrain){
  # varfunc='constant';param.init=NULL;constrain=F
  correlation='compound'
  corstr_ <- ifelse(correlation == 'compound', 2, 
                    ifelse(correlation == 'AR',3,1))
  ## B-splines basis ordinate (x)
  tl <- sl <- 0
  tr <- sr <- 1
  id.unique <- unique(id)
  n <- length(id.unique)
  nbeta <- ndx+deg+1
  nbeta0 <- ndx0+deg0+1
  nbeta1 <- ndx1+deg1+1
  p <- ncol(X)
  nsurv <- ncol(CondProb)
  ## penalty stuff
  p.fac1 <- c(rep(0,deg1), rep(1,ndx1))
  p.fac0 <- c(rep(0,deg0+1), rep(1,ndx0))
  p.fac <- c(rep(0,deg+1), rep(1,ndx))
  lambda=1e1
  Penalty <- diag(c(kronecker(p.fac0,p.fac0),
                    kronecker(p.fac1,kronecker(p.fac,p.fac),FUN='+')*100,
                    p.fac0*1e-8,kronecker(p.fac1,p.fac,FUN='+')*1e-8)) * lambda
  
  
  idx.sts <- which(delta == 1)
  X.sts <- X[idx.sts,]; Y.sts <- Y[idx.sts]; id.sts <- id[idx.sts]
  id.sts.unique <- unique(id.sts)
  idx.lts <- which(delta == 0 & surv == 1)
  X.lts <- X[idx.lts,]; Y.lts <- Y[idx.lts]; id.lts <- id[idx.lts]
  id.lts.unique <- unique(id.lts)
  idx.cen <- which(delta == 0 & surv < 1)
  X.cen <- X[idx.cen,]; Y.cen <- Y[idx.cen]; id.cen <- id[idx.cen]
  id.cen.unique <- unique(id.cen)
  censor <- (length(id.cen.unique)>0) * censor
  if(family == "normal"){
    Bmu01 <- B0[idx.sts,1:nbeta0^2]
    Bmu02 <- B0[idx.lts,1:nbeta0 + nbeta0^2]
    Bmu01.cen <- B0[idx.cen,1:nbeta0^2]
    Bmu02.cen <- B0[idx.cen,1:nbeta0 + nbeta0^2]
    Bmu1 <- matrix(B[idx.sts,1:nbeta^2],ncol=nbeta^2)
    Bmu2 <- matrix(B[idx.lts,1:nbeta + nbeta^2],ncol=nbeta)
    Bmu1.cen <- matrix(B[idx.cen,1:nbeta^2],ncol=nbeta^2)
    Bmu2.cen <- matrix(B[idx.cen,1:nbeta + nbeta^2],ncol=nbeta)
  }else if(family %in% c('poisson','binomial')){
    Bmu01 <- B0[[1]][idx.sts,1:nbeta0^2]
    Bmu02 <- B0[[1]][idx.lts,1:nbeta0 + nbeta0^2]
    Bmu01.cen <- B0[[2]][[1]]
    Bmu02.cen <- B0[[2]][[2]]
    Bmu1 <- matrix(B[[1]][idx.sts,1:nbeta^2],ncol=nbeta^2)
    Bmu2 <- matrix(B[[1]][idx.lts,1:nbeta + nbeta^2],ncol=nbeta)
    Bmu1.cen <- matrix(B[[2]][[1]],ncol=nbeta^2)
    Bmu2.cen <- matrix(B[[2]][[2]],ncol=nbeta)
  }
  if(is.null(param.init)){
    Xs <- X.sts^deg1;Xl <- X.lts^deg1
    if(family == 'normal'){
      theta1.new <- c(coef(lm(Y.sts ~ -1 + Xs + Bmu01)))[1:p]
      theta2.new <- c(coef(lm(Y.lts ~ -1 + Xl + Bmu02)))[1:p]
    }else{
      theta1.new <- c(coef(glm(Y.sts ~ -1 + Xs + Bmu01, family=family)))[1:p]
      theta2.new <- c(coef(glm(Y.lts ~ -1 + Xl + Bmu02, family=family)))[1:p]
    }
  }else{
    theta1.new <- param.init$theta1
    theta2.new <- param.init$theta2
  }
  theta1.new[is.na(theta1.new)] <- 1;theta2.new[is.na(theta2.new)] <- 1
  theta1.new <- theta1.new / sqrt(sum(theta1.new^2))
  theta2.new <- theta2.new / sqrt(sum(theta2.new^2))
  theta.new <- c(theta1.new,theta2.new)
  Xtheta1 <- c(X.sts %*% theta1.new)
  Xtheta2 <- c(X.lts %*% theta2.new)
  
  var.model <- list(A.sts = rep(1, length(idx.sts) + length(idx.lts)),
                    A.lts = rep(1, length(idx.sts) + length(idx.lts)),
                    cen.sts = rep(1, length(idx.cen)),
                    cen.lts = rep(1, length(idx.cen)),
                    rho.sts=0,rho.lts=0)
  
  
  ngamma.sts <- nbeta^2*(nbeta1-1)+nbeta0^2; ngamma.lts <- nbeta*(nbeta1-1)+nbeta0
  gamma1.new <- rep(0,ngamma.sts); gamma2.new <- rep(0,ngamma.lts)
  Dgamma1.cen <- Dgamma2.cen <- 0
  Yhat1.new.test <- rep(0,sum(dat.test$death<=1))
  Yhat2.new.test <- rep(0,sum(dat.test$death>1))
  gamma1.new <- rep(0,ngamma.sts);gamma2.new <- rep(0,ngamma.lts)
  if(is.null(param.init)){
    if(family %in% c('poisson','binomial')){
      B1 <- matrix(Cost_pbase(Xtheta1, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
      B2 <- matrix(Cost_pbase(Xtheta2, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
      Dgamma1 <- cbind(Bmu01, B1[, rep(seq(nbeta1-1), each = nbeta^2)] *
                         Bmu1[, rep(seq(nbeta^2),nbeta1-1)])
      Dgamma2 <- cbind(Bmu02, B2[, rep(seq(nbeta1-1), each = nbeta)] *
                         Bmu2[, rep(seq(nbeta),nbeta1-1)] )
      if(family == 'poisson'){
        Y.sts.tmp <- log(Y.sts+1)
        Y.lts.tmp <- log(Y.lts+1)
      }else if(family=='binomial'){
        Y.sts.tmp <- log((Y.sts+1)/(Y.sts+2))
        Y.lts.tmp <- log((Y.lts+1)/(Y.lts+2))
      }
      gamma1.new <- ginv(t(Dgamma1) %*% Dgamma1) %*% (t(Dgamma1) %*% Y.sts.tmp)
      gamma2.new <- ginv(t(Dgamma2) %*% Dgamma2) %*% (t(Dgamma2) %*% Y.lts.tmp)
    }
  }else{
    gamma1.new <- param.init$gamma1
    gamma2.new <- param.init$gamma2
  }
  gamma.new <- c(gamma1.new,gamma2.new)
  n1 <- 101
  death <- seq(0.01,1.01,length.out = n1)
  nn <- 1:n1#death*100
  dat.test <- data.frame(id = rep(1:n1,nn),
                         surv = rep(death,nn),
                         death = rep(death,nn),
                         time = unlist(sapply(nn,function(x) 1:x)) / 100)
  xt <- matrix(0,nrow(dat.test),p); colnames(xt) <- paste0('x',1:p)
  dat.test <- cbind(dat.test,xt)
  dat.test.sts <- dat.test[which(dat.test$death<=1),]
  dat.test.lts <- dat.test[which(dat.test$death>1),]
  dat.test <- rbind(dat.test.sts,dat.test.lts)
  idx.test.sts <- which(dat.test$death<=1)
  idx.test.lts <- which(dat.test$death>1)
  Bti.test.sts <- Cost_pbase(dat.test.sts$time / dat.test.sts$surv, tl, tr, ndx0, deg0)
  Bsi.test.sts <- Cost_pbase(dat.test.sts$surv, sl, sr, ndx0, deg0)
  Bi.test.sts <- kronecker(matrix(1, ncol = nbeta0), Bti.test.sts) *
    kronecker(Bsi.test.sts, matrix(1, ncol = nbeta0))
  Bi.test.lts <- Cost_pbase(dat.test.lts$time, tl, tr, ndx0, deg0)
  Bi.test <- bdiag(Bi.test.sts, Bi.test.lts)
  # c('comorbidity','raceblack','raceothers,'age75','nochemo')
  X.test.sts <- as.matrix(dat.test[idx.test.sts,4+1:p]) #
  X.test.lts <- as.matrix(dat.test[idx.test.lts,4+1:p]) #
  Bt.test <- Cost_pbase(dat.test$time / dat.test$surv, tl, tr, ndx0, deg0)
  Bs.test <- Cost_pbase(dat.test$surv, sl, sr, ndx0, deg0)
  Bmu01.test <- kronecker(matrix(1, ncol = nbeta0), Bt.test[idx.test.sts,]) *
    kronecker(Bs.test[idx.test.sts,], matrix(1, ncol = nbeta0))
  Bmu02.test <- Bt.test[idx.test.lts,]
  Bmu1.test <- kronecker(matrix(1, ncol = nbeta), Bt.test[idx.test.sts,]) *
    kronecker(Bs.test[idx.test.sts,], matrix(1, ncol = nbeta))
  Bmu2.test <- Bt.test[idx.test.lts,]
  iterate <- 0
  MAX.IT=5
  Penalty <- diag(c(kronecker(p.fac0,p.fac0),
                    kronecker(p.fac1,kronecker(p.fac,p.fac),FUN='+')*100,
                    p.fac0*1e-8,kronecker(p.fac1,p.fac,FUN='+')*1e-8)) * lambda
  repeat{
    gamma <- gamma.new; theta <- theta.new[c(2:p,2:p+p)]
    gamma1 <- gamma1.new; gamma2 <- gamma2.new
    theta1 <- theta1.new; theta2 <- theta2.new
    Yhat1.test <- Yhat1.new.test
    Yhat2.test <- Yhat2.new.test
    
    Xtheta1 <- c(X.sts %*% theta1)
    Xtheta2 <- c(X.lts %*% theta2)
    
    B1 <- matrix(Cost_pbase(Xtheta1, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
    B2 <- matrix(Cost_pbase(Xtheta2, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
    Dgamma1 <- cbind(Bmu01, B1[, rep(seq(nbeta1-1), each = nbeta^2)] *
                       Bmu1[, rep(seq(nbeta^2),nbeta1-1)])
    Dgamma2 <- cbind(Bmu02, B2[, rep(seq(nbeta1-1), each = nbeta)] *
                       Bmu2[, rep(seq(nbeta),nbeta1-1)] )
    y.sts <- Y.sts
    y.lts <- Y.lts
    
    if(family=="poisson"){
      yhat.sts <- c(exp(Dgamma1 %*% gamma1))
      yhat.lts <- c(exp(Dgamma2 %*% gamma2))
      Dgamma1 <- yhat.sts * Dgamma1
      Dgamma2 <- yhat.lts * Dgamma2
      y.sts <- Y.sts - yhat.sts
      y.lts <- Y.lts - yhat.lts
    }else if(family=='binomial'){
      yhat.sts0 <- c(exp(Dgamma1 %*% gamma1))
      yhat.lts0 <- c(exp(Dgamma2 %*% gamma2))
      yhat.sts <- yhat.sts0 / (yhat.sts0+1)
      yhat.lts <- yhat.lts0 / (yhat.lts0+1)
      Dgamma1 <- yhat.sts / (yhat.sts0+1) * Dgamma1
      Dgamma2 <- yhat.lts / (yhat.lts0+1) * Dgamma2
      y.sts <- Y.sts - yhat.sts
      y.lts <- Y.lts - yhat.lts
    }
    Dgamma <- as.matrix(bdiag(Dgamma1,Dgamma2))
    # update gamma
    XtXY.sts <- thetaUpdate( var.model$rho.sts, corstr_, 
                             y.sts, id.sts, id.sts.unique, Dgamma1, 
                             var.model$A.sts, nthread = CORE_NUM )
    XtXY.lts <- thetaUpdate( var.model$rho.lts, corstr_, 
                             y.lts, id.lts, id.lts.unique, Dgamma2, 
                             var.model$A.lts, nthread = CORE_NUM )
    XtX.sts <- XtXY.sts[,1:ngamma.sts]
    XtY.sts <- XtXY.sts[,ngamma.sts+1]
    XtX.lts <- XtXY.lts[,1:ngamma.lts]
    XtY.lts <- XtXY.lts[,ngamma.lts+1]
    XtX <- as.matrix(bdiag(XtX.sts,XtX.lts))
    XtY <- c(XtY.sts,XtY.lts)
    if(censor){
      Xtheta1.cen <- c(X.cen %*% theta1)
      Xtheta2.cen <- c(X.cen %*% theta2)
      B1.cen <- matrix(Cost_pbase(Xtheta1.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
      B2.cen <- matrix(Cost_pbase(Xtheta2.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
      if(family=="normal"){
        Dgamma1.cen <- cbind(Bmu01.cen, B1.cen[, rep(seq(nbeta1-1), each = nbeta^2)] * 
                               Bmu1.cen[, rep(seq(nbeta^2),nbeta1-1)] )
        Dgamma2.cen <- cbind(Bmu02.cen, B2.cen[, rep(seq(nbeta1-1), each = nbeta)] * 
                               Bmu2.cen[, rep(seq(nbeta),nbeta1-1)] )
        # yhat1.cen <- c(Dgamma1.cen %*% gamma1.new)
        # yhat2.cen <- c(Dgamma2.cen %*% gamma2.new)
        y.cen <- Y.cen
      }else if(family=="poisson"){
        ncen <- nrow(B1.cen)
        yhat1.Lst <- matrix(0, nsurv-1, ncen)
        Dgamma1.cen <- matrix(0, ncen, ngamma.sts)
        for(j in 1:(nsurv-1)){
          Bmu01.cen.j <- Bmu01.cen[[j]]
          Bmu1.cen.j <- Bmu1.cen[[j]]
          Dgamma1.cen.j <- cbind(Bmu01.cen.j, B1.cen[, rep(seq(nbeta1-1), each = nbeta^2)] * 
                                   Bmu1.cen.j[, rep(seq(nbeta^2),nbeta1-1)])
          yhat1.cen.j <- CondProb[id.cen, j] * exp(c(Dgamma1.cen.j %*% gamma1))
          yhat1.cen.j[yhat1.cen.j>100] <- 100; yhat1.cen.j[yhat1.cen.j< -100] <- -100
          yhat1.Lst[j,] <- yhat1.cen.j
          Dgamma1.cen <- Dgamma1.cen + yhat1.cen.j * Dgamma1.cen.j
        }
        yhat1.cen <- colSums(yhat1.Lst)
        Dgamma2.cen <- cbind(Bmu02.cen, B2.cen[, rep(seq(nbeta1-1), each = nbeta)] * 
                               Bmu2.cen[, rep(seq(nbeta),nbeta1-1)] )
        yhat2.cen <- CondProb[id.cen, nsurv] * exp(c(Dgamma2.cen %*% gamma2))
        yhat2.cen[yhat2.cen>100] <- 100; yhat2.cen[yhat2.cen< -100] <- -100
        Dgamma2.cen <- yhat2.cen * Dgamma2.cen
        y.cen <- Y.cen - yhat1.cen - yhat2.cen
      }else if(family=="binomial"){
        ncen <- nrow(B1.cen)
        yhat1.Lst <- matrix(0, nsurv-1, ncen)
        yhat1d.Lst <- matrix(0, nsurv-1, ncen)
        Dgamma1.cen <- matrix(0, ncen, ngamma.sts)
        for(j in 1:(nsurv-1)){
          Bmu01.cen.j <- Bmu01.cen[[j]]
          Bmu1.cen.j <- Bmu1.cen[[j]]
          Dgamma1.cen.j <- cbind(Bmu01.cen.j, B1.cen[, rep(seq(nbeta1-1), each = nbeta^2)] * 
                                   Bmu01.cen.j[, rep(seq(nbeta^2),nbeta1-1)])
          yhat1.cen.j0 <- exp(c(Dgamma1.cen.j %*% gamma1))
          yhat1.cen.j <- CondProb[id.cen, j] * yhat1.cen.j0 / (yhat1.cen.j0+1)
          yhat1.cen.j[yhat1.cen.j>100] <- 100; yhat1.cen.j[yhat1.cen.j< -100] <- -100
          yhat1.Lst[j,] <- yhat1.cen.j
          yhat1d.Lst[j,] <- yhat1.cen.j / (yhat1.cen.j0+1)
          Dgamma1.cen <- Dgamma1.cen + yhat1d.Lst[j,] * Dgamma1.cen.j
        }
        yhat1.cen <- colSums(yhat1.Lst)
        Dgamma2.cen <- cbind(Bmu02.cen, B2.cen[, rep(seq(nbeta1-1), each = nbeta)] * 
                               Bmu2.cen[, rep(seq(nbeta),nbeta1-1)] )
        yhat2.cen0 <-  exp(c(Dgamma2.cen %*% gamma2))
        yhat2.cen <- CondProb[id.cen, nsurv] * yhat2.cen0 / (yhat2.cen0+1)
        yhat2.cen[yhat2.cen>100] <- 100; yhat2.cen[yhat2.cen< -100] <- -100
        yhat2d.cen <- yhat2.cen / (yhat2.cen0+1)
        Dgamma2.cen <- yhat2d.cen * Dgamma2.cen
        y.cen <- Y.cen - yhat1.cen - yhat2.cen
      }
      XtXY.cen <- thetaUpdate_cen(var.model$rho.sts, var.model$rho.lts, corstr_, 
                                  y.cen, id.cen, id.cen.unique, 
                                  Dgamma1.cen, Dgamma2.cen, 
                                  var.model$cen.sts, var.model$cen.lts,
                                  nthread = CORE_NUM)
      XtX.cen <- XtXY.cen[, 1:(ngamma.sts+ngamma.lts)]
      XtY.cen <- XtXY.cen[, ngamma.sts+ngamma.lts+1]
    }else{
      XtX.cen <- XtY.cen <- 0
    }
    XtX <- as.matrix((XtX + XtX.cen) / n)
    XtY <- as.matrix((XtY + XtY.cen) / n)
    if(family == "normal"){
      if(!constrain){
        gamma.new <- as.numeric( ginv(XtX + 2*Penalty) %*% (XtY) )
      }else{
        gamma.new <- solve.QP(XtX + 2*Penalty,
                              XtY,
                              diag(c(rep(0,nbeta0^2),
                                     kronecker(rep(1,nbeta^2),c(rep(0,deg1),rep(1,ndx1))),
                                     rep(0,nbeta0),
                                     kronecker(rep(1,nbeta),c(rep(0,deg1),rep(1,ndx1))))),
                              rep(0,ngamma.lts+ngamma.sts))$solution
      }
    }else if(family %in% c('poisson','binomial')){
      gamma.new <- gamma +
        as.numeric( ginv(XtX + 2*Penalty) %*% (XtY - 2* Penalty %*% gamma) )
    }
    
    gamma1.new <- gamma.new[1:ngamma.sts]
    gamma2.new <- gamma.new[1:ngamma.lts+ngamma.sts]
    # update A and rho
    var.model <- VarUpdate(Y.sts, Y.lts, id.sts, id.lts,
                           Dgamma1, Dgamma2, Dgamma1.cen, Dgamma2.cen,
                           gamma1.new, gamma2.new, 
                           deg, ndx, varfunc, corstr_, censor,family=family)
    
    # update theta
    Bp1 <- matrix(Cost_pbase_deriv(Xtheta1, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
    Bp1m <- Bp1[,rep(seq(nbeta1-1), each = nbeta^2)]*
      Bmu1[, rep(seq(nbeta^2),nbeta1-1)] 
    Bp1g <- Bp1m %*% gamma1.new[(nbeta0^2+1):ngamma.sts]
    Jt1 <- rbind(-theta1[2:p] / sqrt(1-sum(theta1[2:p]^2)), diag(rep(1,p-1)))
    Dt1 <- X.sts %*% Jt1 * c(Bp1g)
    Bp2 <- matrix(Cost_pbase_deriv(Xtheta2, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
    Bp2m <- Bp2[, rep(seq(nbeta1-1), each = nbeta)] * 
      Bmu2[, rep(seq(nbeta),nbeta1-1)]
    Bp2g <- Bp2m %*% gamma2.new[(nbeta0+1):ngamma.lts]
    Jt2 <- rbind(-theta2[2:p] / sqrt(1-sum(theta2[2:p]^2)), diag(rep(1,p-1)))
    Dt2 <- X.lts %*% Jt2 * c(Bp2g)
    if(family == "normal"){
      W.sts <- c(Dgamma1 %*% gamma1.new - Dt1 %*% theta1[2:p])
      W.lts <- c(Dgamma2 %*% gamma2.new - Dt2 %*% theta2[2:p])
    }else if(family == 'poisson'){
      Dt1 <- yhat.sts * Dt1
      Dt2 <- yhat.lts * Dt2
      W.sts <- yhat.sts
      W.lts <- yhat.lts
    }else if(family == 'binomial'){
      Dt1 <- yhat.sts / (yhat.sts0 + 1) * Dt1
      Dt2 <- yhat.lts / (yhat.lts0 + 1) * Dt2
      W.sts <- yhat.sts
      W.lts <- yhat.lts
    }
    XtXY.sts <- thetaUpdate( var.model$rho.sts, corstr_, 
                             Y.sts - W.sts, id.sts, id.sts.unique, Dt1, 
                             var.model$A.sts, nthread = CORE_NUM)
    XtXY.lts <- thetaUpdate( var.model$rho.lts, corstr_, 
                             Y.lts - W.lts, id.lts, id.lts.unique, Dt2, 
                             var.model$A.lts, nthread = CORE_NUM)
    XtX.sts <- XtXY.sts[,1:(p-1)]
    XtY.sts <- XtXY.sts[,p]
    XtX.lts <- XtXY.lts[,1:(p-1)]
    XtY.lts <- XtXY.lts[,p]
    XtX <- as.matrix(bdiag(XtX.sts,XtX.lts))
    XtY <- c(XtY.sts,XtY.lts)
    
    
    if(censor){
      if(family=='normal'){
        Bp1.cen <- matrix(Cost_pbase_deriv(Xtheta1.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
        Bp1m <- Bp1.cen[, rep(seq(nbeta1-1), each = nbeta^2)] * 
          Bmu1.cen[, rep(seq(nbeta^2), nbeta1-1)]
        Bp1g <- c(Bp1m %*% gamma1.new[(nbeta0^2+1):ngamma.sts])
        Dt1.cen <- Bp1g * (X.cen %*% Jt1)
        Bp2.cen <- matrix(Cost_pbase_deriv(Xtheta2.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
        Bp2m <- Bp2.cen[, rep(seq(nbeta1-1), each = nbeta)] * 
          Bmu2.cen[, rep(seq(nbeta),nbeta1-1)]
        Bp2g <- c(Bp2m %*% gamma2.new[(nbeta0+1):ngamma.lts])
        Dt2.cen <- Bp2g * (X.cen %*% Jt2)
        W.cen <- c(Dgamma1.cen %*% gamma1.new - Dt1.cen %*% theta1[2:p]) +
          c(Dgamma2.cen %*% gamma2.new - Dt2.cen %*% theta2[2:p])
      }else if(family=="poisson"){
        Bp1.cen <- matrix(Cost_pbase_deriv(Xtheta1.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
        Dt1.cen <- matrix(0, ncen, p-1)
        for(j in 1:(nsurv-1)){
          Bmu01.cen.j <- Bmu01.cen[[j]]
          Bmu1.cen.j <- Bmu1.cen[[j]]
          Bp1m.j <- Bp1.cen[, rep(seq(nbeta1-1), each = nbeta^2)] * 
            Bmu1.cen.j[, rep(seq(nbeta^2),nbeta1-1)]
          Bp1g.j <- c(Bp1m.j %*% gamma1.new[(nbeta0^2+1):ngamma.sts])
          Dt1.cen <- Dt1.cen + yhat1.Lst[j,] * Bp1g.j * (X.cen %*% Jt1)
        }
        Bp2.cen <- matrix(Cost_pbase_deriv(Xtheta2.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
        Bp2m <- Bp2.cen[, rep(seq(nbeta1-1), each = nbeta)] * 
          Bmu2.cen[, rep(seq(nbeta),nbeta1-1)]
        Bp2g <- c(Bp2m %*% gamma2.new[(nbeta0+1):ngamma.lts])
        Dt2.cen <- yhat2.cen * Bp2g * (X.cen %*% Jt2)
        W.cen <- yhat1.cen + yhat2.cen
      }else if(family=="binomial"){
        Bp1.cen <- matrix(Cost_pbase_deriv(Xtheta1.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
        Dt1.cen <- matrix(0, ncen, p-1)
        for(j in 1:(nsurv-1)){
          Bmu1.cen.j <- Bmu1.cen[[j]]
          Bmu01.cen.j <- Bmu01.cen[[j]]
          Bp1m.j <- Bp1.cen[, rep(seq(nbeta1-1), each = nbeta^2)] * 
            Bmu1.cen.j[, rep(seq(nbeta^2),nbeta1-1)]
          Bp1g.j <- c(Bp1m.j %*% gamma1.new[(nbeta0^2+1):ngamma.sts])
          Dt1.cen <- Dt1.cen + yhat1d.Lst[j,] * Bp1g.j * (X.cen %*% Jt1)
        }
        Bp2.cen <- matrix(Cost_pbase_deriv(Xtheta2.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
        Bp2m <- Bp2.cen[, rep(seq(nbeta1-1), each = nbeta)] * 
          Bmu2.cen[, rep(seq(nbeta),nbeta1-1)]
        Bp2g <- c(Bp2m %*% gamma2.new[(nbeta0+1):ngamma.lts])
        Dt2.cen <- yhat2d.cen * Bp2g * (X.cen %*% Jt2)
        W.cen <- yhat1.cen + yhat2.cen
      }
      W.cen[W.cen>50] <- 50; W.cen[W.cen< -50] <- -50
      XtXY.cen <- thetaUpdate_cen(var.model$rho.sts, var.model$rho.lts, corstr_, 
                                  Y.cen - W.cen, id.cen, id.cen.unique, 
                                  Dt1.cen, Dt2.cen, 
                                  var.model$cen.sts, var.model$cen.lts,
                                  nthread = CORE_NUM)
      
      XtX.cen <- XtXY.cen[, 1:(2*(p-1))]
      XtY.cen <- XtXY.cen[, 2*p-1]
    }else{
      XtX.cen <- XtY.cen <- 0
    }
    XtX <- as.matrix((XtX + XtX.cen) / n)
    XtY <- as.matrix((XtY + XtY.cen) / n)
    if(family == "normal"){
      theta.new <- as.numeric( ginv(XtX) %*% (XtY) )
    }else if(family  %in% c('poisson','binomial')){
      theta.new <- theta + 
        as.numeric( ginv(XtX) %*% (XtY) )
    }
    
    theta.new[(theta.new>1) | (theta.new < -1)] <- 0
    theta1.new <- theta.new[1:(p-1)];theta2.new <- theta.new[1:(p-1)+(p-1)]
    st1 <- sum(theta1.new^2)
    st2 <- sum(theta2.new^2)
    if(st1>1) theta1.new <- theta1.new / sqrt(st1+1e-6)
    if(st2>1) theta2.new <- theta2.new / sqrt(st2+1e-2)
    theta1.new <- c(sqrt(1-sum(theta1.new^2)),theta1.new)
    theta2.new <- c(sqrt(1-sum(theta2.new^2)),theta2.new)
    theta.new <- c(theta1.new, theta2.new)
    ########################plot
    
    
    Xtheta1.test <- c(X.test.sts %*% theta1.new)
    Xtheta2.test <- c(X.test.lts %*% theta2.new)
    B1.test <- matrix(Cost_pbase(Xtheta1.test, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
    B2.test <- matrix(Cost_pbase(Xtheta2.test, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
    Dgamma1.test <- cbind(Bmu01.test, B1.test[, rep(seq(nbeta1-1), each = nbeta^2)] * 
                            Bmu1.test[, rep(seq(nbeta^2),nbeta1-1)])
    Dgamma2.test <- cbind(Bmu02.test, B2.test[, rep(seq(nbeta1-1), each = nbeta)] * 
                            Bmu2.test[, rep(seq(nbeta),nbeta1-1)] )
    
    Yhat1.new.test <- Dgamma1.test %*% gamma1.new
    Yhat2.new.test <- Dgamma2.test %*% gamma2.new
    par(mfrow = c(1,2))
    plot(NULL, xlim = c(0,1),ylim = c(0,.2),xlab='Time', ylab='Y', 
         main=paste(c('Iteration #',iterate,'\n theta1=(',round(theta1.new,2),')\n theta2=(',round(theta2.new,2),')'),collapse =' '))
    
    for(dd in c(0.4,0.6,0.8,0.9,1)){
      idx <- which(dat.test.sts$death==dd)
      lines(Yhat1.new.test[idx]~dat.test.sts$time[idx],col=1,lty=1)
    }
    lines(Yhat2.new.test~dat.test.lts$time,col=1,lty=2)
    
    plot(c(gamma.new,theta.new[c(2:p,2:p+p)])-
           c(gamma,theta) , ylab='difference of f(x)',
         main = paste0('Iteration #',iterate))
    ##########################end plot
    iterate <- iterate + 1
    if( Convergence(c(gamma.new,theta.new[c(2:p,2:p+p)]),
                    c(gamma,theta)) | iterate > MAX.IT){
      converge <- 1
      break
    }
  };iterate=0
  
  
   return( list(
    gamma1 = gamma1, gamma2 = gamma2,
    theta1 = theta1, theta2 = theta2, varmodel = var.model) )
}


# 
# if(is.null(theta.old)){
#  lambda0 <- lambda
#  if(is.null(lambda)) lambda0 <- 1e-6
#  Penalty <- lambda0 * diag(p.fac)
#  
#  theta.old <- as.matrix(ginv( (t(B) %*% B) / n + 2 * Penalty) %*% 
#              ( (t(B) %*% Y) / n ))
#  
#  plot(NULL, xlim = c(0,1), ylim = c(0,20),
#    main = "Initial theta", xlab='time', ylab='Y')
#  for(i in seq(.2,1,.2)) lines(dat.test$time[which(dat.test$death==i)],
#                Bi[which(dat.test$death==i),] %*% 
#                 theta.old[c(1:nbeta^2, p*nbeta^2+1:nbeta)], 
#                col = 1+i*10,lty=1) # STS
#  lines(dat.test$time[which(dat.test$death>1)], 
#     Bi[which(dat.test$death>1),] %*% 
#      theta.old[c(1:nbeta^2, p*nbeta^2+1:nbeta)], col=2,lty=2) # LTS
# }
# theta.init <- theta.old
# ntheta <- length(theta.init)
# idx.theta.sts <- 1:(p*nbeta^2)
# idx.theta.lts <- 1:(p*nbeta)+(p*nbeta^2)
# ntheta.sts <- length(idx.theta.sts)
# ntheta.lts <- length(idx.theta.lts)
# 
# idx.sts <- which(delta == 1)
# idx.lts <- which(delta == 0 & s == 1)
# idx.cen <- which(delta == 0 & s < 1)
# 
# id.sts <- id[idx.sts]
# id.lts <- id[idx.lts]
# id.cen <- id[idx.cen]
# 
# id.sts.unique <- unique(id.sts)
# id.lts.unique <- unique(id.lts)
# id.cen.unique <- unique(id.cen)
# 
# ## penalty stuff
# 
# if(is.null(lambda)){
#  lambdas <- 10^c(-4,4)
#  lambda <- select_cost_lambda(Bh, K, eo, Ytilde, w, DtD,
#                beta.new = NULL, 
#                MON, TOL1, MAX.IT, ndx, deg, n,
#                lambdas)
# }
# 
# Penalty <- lambda * diag(p.fac)
# Penalty_var <- lambda * diag(c(rep(0,deg+1),rep(1,ndx)))
# 
# theta <- fit_cost(B, theta.old, z, id, varfunc, 
#          idx.sts, idx.lts,
#          id.sts, id.lts,
#          id.sts.unique, id.lts.unique, id.cen.unique,
#          idx.theta.sts, idx.theta.lts,
#          Penalty, Penalty_var, Bi, dat.test, correlation,
#          MON,TOL1, MAX.IT, ndx, deg, n,
#          ntheta.sts, ntheta.lts, ntheta)

VarUpdate <- function(Y.sts, Y.lts, id.sts, id.lts,
                      Dgamma1, Dgamma2, Dgamma1.cen, Dgamma2.cen,
                      gamma1, gamma2, 
                      deg, ndx, varfunc, corstr_, censor,family='normal'){
  var.sts=var.lts=cen.sts=cen.lts=NULL
  CORE_NUM <- parallel::detectCores()-5
  p.fac <- c(rep(0,deg+1),rep(1,ndx))
  Penalty_var <- diag(p.fac) * 0.01
  
  Yhat.sts <- c(Dgamma1 %*% gamma1)
  Yhat.lts <- c(Dgamma2 %*% gamma2)
  if(family == 'poisson'){
    Yhat.sts <- exp(Yhat.sts)
    Yhat.lts <- exp(Yhat.lts)
  }else if(family == 'binomial'){
    Yhat.sts <- exp(Yhat.sts) / (exp(Yhat.sts)+1)
    Yhat.lts <- exp(Yhat.lts) / (exp(Yhat.lts)+1)
  }
  Yhat.sts[Yhat.sts>100] <- 100
  Yhat.lts[Yhat.lts>100] <- 100
  id.sts.unique <- unique(id.sts)
  id.lts.unique <- unique(id.lts)
  n.sts <- length(id.sts)
  n.lts <- length(id.lts)
  
  err.sts <- Y.sts - Yhat.sts
  err.lts <- Y.lts - Yhat.lts
  
  logerr.sq.sts = log(err.sts^2)
  logerr.sq.lts = log(err.lts^2)
  
  if(varfunc == 'constant'){
    A.sts <- rep(1,length(Yhat.sts))
    A.lts <- rep(1,length(Yhat.lts))
    if(censor){
      cen.sts <- rep(1,nrow(Dgamma1.cen))
      cen.lts <- rep(1,nrow(Dgamma2.cen))
    }
  }else if(varfunc == 'linear'){
    varmodel.sts <- lm(logerr.sq.sts ~ Yhat.sts)
    varmodel.lts <- lm(logerr.sq.lts ~ Yhat.lts)
    A.sts <- sqrt(c(exp(predict(varmodel.sts, 
                                newdata=data.frame(Yhat.sts=Yhat.sts)))))
    A.lts <- sqrt(c(exp(predict(varmodel.lts, 
                                newdata=data.frame(Yhat.lts=Yhat.lts)))))
    if(censor){
      Yhat1 <- c(Dgamma1.cen %*% gamma1)
      Yhat2 <- c(Dgamma2.cen %*% gamma2)
      if(family == 'poisson'){
        Yhat1 <- exp(Yhat1)
        Yhat2 <- exp(Yhat2)
      }else if(family == 'binomial'){
        Yhat1 <- exp(Yhat1) / (exp(Yhat1)+1)
        Yhat2 <- exp(Yhat2) / (exp(Yhat2)+1)
        Yhat1[Yhat1>1] <- 1
        Yhat2[Yhat2>1] <- 1
      }
      
      Yhat1[Yhat1>100] <- 100
      Yhat2[Yhat2>100] <- 100
      cen.sts <- sqrt(c(exp(predict(varmodel.sts, 
                                    newdata=data.frame(Yhat.sts=Yhat1)))))
      cen.lts <- sqrt(c(exp(predict(varmodel.lts, 
                                    newdata=data.frame(Yhat.lts=Yhat2)))))
    }
  }else{
    B.Yhat.sts <- Cost_pbase(Yhat.sts, xl=0, xr=20, ndx, deg)
    B.Yhat.lts <- Cost_pbase(Yhat.lts, xl=0, xr=20, ndx, deg)
    
    var.sts <- as.matrix(ginv( (t(B.Yhat.sts) %*% B.Yhat.sts) / n.sts + 2 * Penalty_var) %*% 
                           ( (t(B.Yhat.sts) %*% logerr.sq.sts) / n.sts ))
    var.lts <- as.matrix(ginv( (t(B.Yhat.lts) %*% B.Yhat.lts) / n.lts + 2 * Penalty_var) %*% 
                           ( (t(B.Yhat.lts) %*% logerr.sq.lts) / n.lts ))
    A.sts <- sqrt(exp(B.Yhat.sts %*% var.sts))
    A.lts <- sqrt(exp(B.Yhat.lts %*% var.lts))
    if(censor){
      Yhat1 <- c(Dgamma1.cen %*% gamma1)
      Yhat2 <- c(Dgamma2.cen %*% gamma2)
      if(family == 'poisson'){
        Yhat1 <- exp(Yhat1)
        Yhat2 <- exp(Yhat2)
      }else if(family == 'binomial'){
        Yhat1 <- exp(Yhat1) / (exp(Yhat1)+1)
        Yhat2 <- exp(Yhat2) / (exp(Yhat2)+1)
        Yhat1[Yhat1>1] <- 1
        Yhat2[Yhat2>1] <- 1
      }
      
      Yhat1[Yhat1>100] <- 100
      Yhat2[Yhat2>100] <- 100
      Bcen1.Yhat <- Cost_pbase(Yhat1, xl=0, xr=20, ndx, deg)
      Bcen2.Yhat <- Cost_pbase(Yhat2, xl=0, xr=20, ndx, deg)
      cen.sts <- sqrt(exp(Bcen1.Yhat %*% var.sts))
      cen.lts <- sqrt(exp(Bcen2.Yhat %*% var.lts))
    }
  }
  rho.sts <- thetaCorr(err.sts, id.sts, id.sts.unique, corstr_, nthread = CORE_NUM)
  rho.lts <- thetaCorr(err.lts, id.lts, id.lts.unique, corstr_, nthread = CORE_NUM)
  return(list(A.sts=A.sts, A.lts=A.lts,
              var.sts=var.sts, var.lts=var.lts,
              cen.sts=cen.sts, cen.lts=cen.lts,
              rho.sts=rho.sts, rho.lts=rho.lts))
}