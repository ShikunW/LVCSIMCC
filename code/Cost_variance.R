library(mgcv)
Cost_variance <- function(Y, X, B0, B, id, delta, surv, surv.out,
                          ndx, deg,  ndx0, deg0, ndx1, deg1, param,
                          lambda, MON, correlation, 
                          varfunc='constant', censor=F, family = 'normal'){
  # varfunc='spline';
  gamma1 <- param$gamma1; gamma2 <- param$gamma2
  theta1 <- param$theta1; theta2 <- param$theta2
  CondProb <- surv.out$cond_surv$prob
  correlation='compound'
  corstr_ <- ifelse(correlation == 'compound', 2, ifelse(correlation == 'AR',3,1))
  CORE_NUM <- parallel::detectCores()-5
  ## B-splines basis ordinate (x)
  tl <- sl <- 0
  tr <- sr <- 1
  nbeta <- ndx+deg+1
  nbeta0 <- ndx0+deg0+1
  nbeta1 <- ndx1+deg1+1
  p <- ncol(X)
  nparam <- nbeta0-2+p
  nsurv <- ncol(CondProb)
  ## penalty stuff
  
  # p.fac1 <- c(rep(0,deg1), rep(1,ndx1))
  # p.fac0 <- c(rep(0,deg0+1), rep(1,ndx0))
  # p.fac <- c(rep(0,deg+1), rep(1,ndx))
  # Penalty <- diag(c(kronecker(p.fac0,p.fac0),p.fac0,p.fac1,p.fac1)) * lambda
  

  idx.sts <- which(delta == 1)
  X.sts <- X[idx.sts,]; Y.sts <- Y[idx.sts]; id.sts <- id[idx.sts]
  id.sts.unique <- unique(id.sts)
  idx.lts <- which(delta == 0 & surv == 1)
  X.lts <- X[idx.lts,]; Y.lts <- Y[idx.lts]; id.lts <- id[idx.lts]
  id.lts.unique <- unique(id.lts)
  idx.cen <- which(delta == 0 & surv < 1)
  X.cen <- X[idx.cen,]; Y.cen <- Y[idx.cen]; id.cen <- id[idx.cen]
  Y <- c(Y.sts,Y.lts,Y.cen)
  id.cen.unique <- unique(id.cen)
  
  censor <- (length(id.cen.unique)>0) * censor
  if(!censor){
    Y <- c(Y.sts, Y.lts)
    id <- c(id.sts, id.lts)
  }
  id.unique <- unique(id)
  n <- length(id.unique)
  ngamma.sts <- nbeta^2*(nbeta1-1)+nbeta0^2; ngamma.lts <- nbeta*(nbeta1-1)+nbeta0
  ntheta <- ngamma.sts + ngamma.lts
  
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
  
  # gamma <- gamma; theta <- theta
  # gamma1 <- gamma1; gamma2 <- gamma2
  # theta1 <- theta1; theta2 <- theta2
  
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
  if(family=="normal"){
    yhat.sts <- c(Dgamma1 %*% gamma1)
    yhat.lts <- c(Dgamma2 %*% gamma2)
  }else if(family=="poisson"){
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
  
  Dgamma1.cen <- Dgamma2.cen <- Dt1.cen <- Dt2.cen <- NULL
  Dbeta.cen <- matrix(0,ntheta,nparam)
  
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
      yhat1.cen <- c(Dgamma1.cen %*% gamma1)
      yhat2.cen <- c(Dgamma2.cen %*% gamma2)
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
      
      Dbeta.cen <- matrix(0,ntheta,nparam)
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
      
      Dbeta.cen <- matrix(0,ntheta,nparam)
    }
  }
  
  Bp1 <- matrix(Cost_pbase_deriv(Xtheta1, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
  Bp1m <- Bp1[,rep(seq(nbeta1-1), each = nbeta^2)]*
    Bmu1[, rep(seq(nbeta^2),nbeta1-1)] 
  Bp1g <- Bp1m %*% gamma1[(nbeta0^2+1):ngamma.sts]
  Jt1 <- rbind(-theta1[2:p] / sqrt(1-sum(theta1[2:p]^2)), diag(rep(1,p-1)))
  Dt1 <- X.sts %*% Jt1 * c(Bp1g)
  Bp2 <- matrix(Cost_pbase_deriv(Xtheta2, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
  Bp2m <- Bp2[, rep(seq(nbeta1-1), each = nbeta)] * 
    Bmu2[, rep(seq(nbeta),nbeta1-1)]
  Bp2g <- Bp2m %*% gamma2[(nbeta0+1):ngamma.lts]
  Jt2 <- rbind(-theta2[2:p] / sqrt(1-sum(theta2[2:p]^2)), diag(rep(1,p-1)))
  Dt2 <- X.lts %*% Jt2 * c(Bp2g)
  
  if(family=="poisson"){
    Dt1 <- yhat.sts * Dt1
    Dt2 <- yhat.lts * Dt2
  }else if(family == 'binomial'){
    Dt1 <- yhat.sts / (yhat.sts0 + 1) * Dt1
    Dt2 <- yhat.lts / (yhat.lts0 + 1) * Dt2
  }
  
  if(censor){
    if(family=='normal'){
      Bp1.cen <- matrix(Cost_pbase_deriv(Xtheta1.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
      Bp1m <- Bp1.cen[, rep(seq(nbeta1-1), each = nbeta^2)] * 
        Bmu1.cen[, rep(seq(nbeta^2), nbeta1-1)]
      Bp1g <- Bp1m %*% gamma1[(nbeta0^2+1):ngamma.sts]
      Dt1.cen <- X.cen %*% Jt1 * c(Bp1g)
      Bp2.cen <- matrix(Cost_pbase_deriv(Xtheta2.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
      Bp2m <- Bp2.cen[, rep(seq(nbeta1-1), each = nbeta)] * 
        Bmu2.cen[, rep(seq(nbeta),nbeta1-1)]
      Bp2g <- Bp2m %*% gamma2[(nbeta0+1):ngamma.lts]
      Dt2.cen <- X.cen %*% Jt2 * c(Bp2g)
    }else if(family=="poisson"){
      Bp1.cen <- matrix(Cost_pbase_deriv(Xtheta1.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
      Dt1.cen <- matrix(0, ncen, p-1)
      for(j in 1:(nsurv-1)){
        Bmu01.cen.j <- Bmu01.cen[[j]]
        Bmu1.cen.j <- Bmu1.cen[[j]]
        Bp1m.j <- Bp1.cen[, rep(seq(nbeta1-1), each = nbeta^2)] * 
          Bmu1.cen.j[, rep(seq(nbeta^2),nbeta1-1)]
        Bp1g.j <- c(Bp1m.j %*% gamma1[(nbeta0^2+1):ngamma.sts])
        Dt1.cen <- Dt1.cen + yhat1.Lst[j,] * Bp1g.j * (X.cen %*% Jt1)
      }
      Bp2.cen <- matrix(Cost_pbase_deriv(Xtheta2.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
      Bp2m <- Bp2.cen[, rep(seq(nbeta1-1), each = nbeta)] * 
        Bmu2.cen[, rep(seq(nbeta),nbeta1-1)]
      Bp2g <- c(Bp2m %*% gamma2[(nbeta0+1):ngamma.lts])
      Dt2.cen <- yhat2.cen * Bp2g * (X.cen %*% Jt2)
    }else if(family=="binomial"){
      Bp1.cen <- matrix(Cost_pbase_deriv(Xtheta1.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
      Dt1.cen <- matrix(0, ncen, p-1)
      for(j in 1:(nsurv-1)){
        Bmu1.cen.j <- Bmu1.cen[[j]]
        Bmu01.cen.j <- Bmu01.cen[[j]]
        Bp1m.j <- Bp1.cen[, rep(seq(nbeta1-1), each = nbeta^2)] * 
          Bmu1.cen.j[, rep(seq(nbeta^2),nbeta1-1)]
        Bp1g.j <- c(Bp1m.j %*% gamma1[(nbeta0^2+1):ngamma.sts])
        Dt1.cen <- Dt1.cen + yhat1d.Lst[j,] * Bp1g.j * (X.cen %*% Jt1)
      }
      Bp2.cen <- matrix(Cost_pbase_deriv(Xtheta2.cen, sl, sr, ndx1, deg1)[,-1], ncol=nbeta1-1)
      Bp2m <- Bp2.cen[, rep(seq(nbeta1-1), each = nbeta)] * 
        Bmu2.cen[, rep(seq(nbeta),nbeta1-1)]
      Bp2g <- c(Bp2m %*% gamma2[(nbeta0+1):ngamma.lts])
      Dt2.cen <- yhat2d.cen * Bp2g * (X.cen %*% Jt2)
    }
  }
  # update A and rho
  var.model <- VarUpdate(Y.sts, Y.lts, id.sts, id.lts,
                         Dgamma1, Dgamma2, Dgamma1.cen, Dgamma2.cen, 
                         gamma1, gamma2, 
                         deg, ndx, varfunc, corstr_, censor)
  
  Yhat.sts <- yhat.sts
  Yhat.lts <- yhat.lts
  Yhat.cen <- NULL
  if(censor){
    Yhat.cen <- yhat1.cen + yhat2.cen
  }
  
  Yhat <- c(Yhat.sts, Yhat.lts, Yhat.cen)
  err <- Y - Yhat
  
  D1 <- cbind(Dgamma1, Dt1)
  D2 <- cbind(Dgamma2, Dt2)
  D1.cen <- cbind(Dgamma1.cen, Dt1.cen)
  D2.cen <- cbind(Dgamma2.cen, Dt2.cen)
  # D1 <- cbind(Dt1)
  # D2 <- cbind(Dt2)
  # D1.cen <- cbind(Dt1.cen)
  # D2.cen <- cbind(Dt2.cen)
  DD1 <- rbind(D1,matrix(0,nrow(D2),ncol(D1)),D1.cen)
  DD2 <- rbind(matrix(0,nrow(D1),ncol(D2)),D2,D2.cen)
  A.sts <- c(var.model$A.sts, rep(0, length(var.model$A.lts)), var.model$cen.sts)
  A.lts <- c(rep(0, length(var.model$A.sts)), var.model$A.lts, var.model$cen.lts)
  
  ntheta <- ncol(DD1) + ncol(DD2)
  ntheta0 <- ncol(Dgamma1) + ncol(Dgamma2)
  
  
  
  gradn_beta <- matrix(0,nparam,n) # surv.out$var_surv$gradn # 
  deriv_beta <- matrix(0,nrow(DD1),nparam)
  XtXY <- thetaVariance(var.model$rho.sts, var.model$rho.lts, corstr_,
                        id, id.unique, DD1, DD2, gradn_beta, deriv_beta, 
                        A.sts, A.lts, 
                        err, ntheta, nparam, nthread = CORE_NUM)   
  
  bread <- meat <- matrix(0, ntheta+nparam, ntheta+nparam)
  bread[1:ntheta,] <- XtXY[1:ntheta,]
  bread[ntheta+1:nparam, ntheta+1:nparam] <- surv.out$var_surv$hess
  
  Pmat <- as.matrix(bdiag(Penalty[1:ngamma.sts,1:ngamma.sts],matrix(0,p-1,p-1),
                          Penalty[ngamma.sts+1:ngamma.lts,ngamma.sts+1:ngamma.lts],matrix(0,p-1,p-1)))
  # bread[1:ntheta,1:ntheta] <- bread[1:ntheta,1:ntheta] + 2 * Pmat * lambda
  bread[1:ntheta,1:ntheta] <- bread[1:ntheta,1:ntheta] + 2 * Pmat
  # bread[1:nbeta + ntheta, 1:nbeta + ntheta] <- surv.out$var_surv$hess
  meat <- XtXY[ntheta+1:(ntheta+nparam),]
  # bread[bread>1e5] <- 1; meat[meat>1e5] <- 1
  bread_inv <- solve(bread + 1e-6*diag(ntheta+nparam))#+ 1e-6*diag(ntheta)
  
  variance <- (bread_inv %*% meat %*% bread_inv)[1:ntheta,1:ntheta]
  g1 <- diag(nbeta0^2+nbeta^2*(nbeta1-1))
  g2 <- cbind(-theta1[2:p] / sqrt(1-sum(theta1[2:p]^2)), diag(p-1))
  g3 <- diag(nbeta0+nbeta*(nbeta1-1))
  g4 <- cbind(-theta2[2:p] / sqrt(1-sum(theta2[2:p]^2)), diag(p-1))
  g <- bdiag(g1,g2,g3,g4)
  v1 <- t(g) %*% variance %*% g
  # plot(diag(variance))
  var.model=list(rho.sts=var.model$rho.sts,rho.lts=var.model$rho.lts)
  return(list(variance = v1,varmodel = var.model))
}
