# surv0 <- runif(100)
# delta0 <- rbinom(100,1,.5)
# X0 <- matrix(runif(200),ncol=2)

library(MASS)
library(Matrix)
Cost_survival0 <- function(surv0, delta0, X0, id0,
                           ndx0 = 4, deg0 = 2, lambda = 1e-5){
  MON=T; TOL1=1e-5; MAX.IT=10
  ntau <- 101
  ndx <- ndx0; deg <- deg0
  tau <- seq(0,1,length.out = ntau)
  Y <- surv0
  n <- length(delta0)
  nbeta0 <- ndx+deg+1
  p <- ncol(X0)
  nparam <- nbeta0+p
  K <- apply(outer(Y, tau, FUN = function(x,y) 1*(x>y)),1,sum) + 1
  Ytilde <- lapply(1:n, function(i) c(rep(0,K[i]-1), delta0[i]))
  # o <- lapply(1:n, function(i) c(min(tau[1], Y[i]), )
  
  o <- outer(Y, c(tau,2), FUN=function(x,y) pmin(x,y))
  o1 <- o[,1:(ntau-1)]
  o2 <- o[,3:(ntau+1)]
  o <- apply(log((o2-o1)/2),1,function(x) x[x>-30])
  eo <- lapply(1:n, function(i) c(pmin(tau[2],Y[i]),exp(o[[i]]))) # exp(o_ik)
  
  Bh <- Cost_bbase(tau, 0, 1, ndx, deg)
  # lambda <- 0
  D <- diag(ndx + deg+1)
  DtD <- diag(c(rep(0,deg),rep(1,ndx+1),rep(0,p)))
  
  n2 <- 101; surv00 <- seq(0,1,length.out = n2)
  w <- Cost_bbase(surv00, 0, 1, ndx, deg)
  
  if(is.null(lambda)){
    lambdas <- 10^c(-6:6)
    lambda <- select_survival_lambda(Bh, K, eo, Ytilde, w, DtD,
                                     beta.new = NULL, 
                                     MON, TOL1, MAX.IT, ndx, deg, n,
                                     lambdas)
  }
  
  P <- 2*lambda * DtD
  beta <- fit_survival(Bh, K, eo, Ytilde, w, P, X0,
                       beta.new = NULL,ndx, deg, p, n)
  
  marg_surv <- marginal_survival(tau, Bh, eo,  
                                 surv0, delta0, K,
                                 beta, ndx, deg, p)
  
  cond_surv <- conditional_survival(marg_surv, n, nparam, nbeta0, X0, id0,
                                    surv0, delta0)
  
  var_surv <- variance_survival(n, K, Bh, Ytilde, eo, P, beta,
                                surv0, delta0, ndx, deg, p)
  # beta_var_meat <- Reduce('+', lapply(1:n,function(i) tcrossprod(var_surv$gradn[,i])))
  # beta_var_bread <- var_surv$hess
  # beta_var <- t(beta_var_bread) %*% beta_var_meat %*% beta_var_bread 
  
  return(list(beta=beta, surv=cond_surv$surv,prob=cond_surv$prob, var_surv=var_surv))
}


select_survival_lambda <- function(Bh, K, eo, Ytilde, w, DtD,
                                   beta.new = NULL, 
                                   ndx, deg, n,
                                   lambdas = c(0,10^seq(-6,6,1))){
  params <- c(); Devs <- c(); 
  n2 <- 101; surv00 <- seq(0,1,length.out = n2)
  TOL1 <- 1e-5
  
  par(mfrow=c(1,2))
  if(MON){
    plot(NULL, xlim = c(0,1), ylim = c(-6,2), 
         xlab = 'Time', ylab = 'log-hazard', 
         main = 'Estimated log-hazard function')
    logHazs <- c()
  }
  for(lambda in lambdas){
    P <- lambda * DtD
    
    beta <- fit_survival(Bh, K, eo, Ytilde, w, P,
                         beta.new, MON=F, TOL1, MAX.IT, ndx, deg, n)
    
    Bhb <- Bh %*% beta
    lmu <- lapply(1:n, function(i) Bhb[1:K[i]] + log(eo[[i]]))
    Dev <- Reduce('+', lapply(1:n, function(i) sum(Ytilde[[i]] * lmu[[i]] + lmu[[i]])))
    
    params <- cbind(params, param)
    Devs <- c(Devs, Dev)
    if(MON){
      par(new=TRUE)
      logHaz <- w %*% beta
      plot((logHaz) ~ surv00, xlab = '', ylab = '', 
           xlim = c(0,1), ylim = c(-6,2), type='l')
      logHazs <- cbind(logHazs, logHaz)
    }
  }
  s <- sqrt(apply(params^2,2,sum)); s <- s / s[1]
  p <- ndx + deg
  n <- length(surv0)
  gcv <- (Devs / (n * (1-p*s/n)^2))
  lambda.choose <- lambdas[2:length(lambdas)][which.min(gcv)]
  
  if(MON){
    par(new=TRUE)
    logHaz <- logHazs[,2:length(lambdas)][,which.min(gcv)]
    plot((logHaz) ~ surv00, xlab = '', ylab = '', 
         xlim = c(0,1), ylim = c(-6,2), type='l', col=2)
    plot(gcv~log10(lambdas), type = 'l')
    abline(v = log10(lambda.choose),col=2)
  }
  
  return(lambda.choose)
}


fit_survival <- function(Bh, K, eo, Ytilde, w, P, X0,
                         beta.new = NULL,ndx, deg, p, n){
  TOL1 <- 1e-5
  if(is.null(beta.new)) beta.new <- rep(1,ndx+deg+p+1)
  n2 <- 101; surv00 <- seq(0,1,length.out = n2)
  MON=T; TOL1=1e-5; MAX.IT=10
  tol <- 1; niter <- 0
  while(tol > TOL1 && niter < MAX.IT){
    beta.old <- beta.new
    beta0.old <- beta.old[1:(ndx+deg+1)]
    beta1.old <- beta.old[ndx+deg+1 + 1:p]
    eBhb <- exp(Bh %*% beta0.old)
    eXb <- exp(X0 %*% beta1.old)
    rate <- lapply(1:n, function(i) eBhb[1:K[i]] * eo[[i]] * eXb[i])
    hess <-  -Reduce('+', lapply(1:n, function(i) 
      Reduce('+',lapply(1:K[i],function(j) 
        tcrossprod(c(Bh[j,],X0[i,])) * rate[[i]][j]))) ) 
    gradn <-  Reduce('+', lapply(1:n, function(i) 
      Reduce('+',lapply(1:K[i],function(j) 
        c(Bh[j,],X0[i,]) * (Ytilde[[i]][j] - rate[[i]][j])))) ) 
    beta.new <- beta.old - ginv(hess-P) %*% (gradn - P %*% beta.old) # Newton-Raphson
    tol <- max(abs(beta.old - beta.new) / abs(beta.old + 1e-10))
    niter <- niter + 1
  }
  
  if(MON){
    logHaz <- w %*% beta.new[1:(ndx+deg+1)]
    plot((logHaz) ~ surv00, type='l', 
         xlab = 'Time', ylab = 'log-hazard', 
         main = 'Estimated log-hazard function')
  }
  
  if(niter > (MAX.IT-1)) {
    warning(paste("parameter estimates did NOT converge in",
                  MAX.IT,
                  "iterations. Increase MAX.IT in control."))
  }
  return(beta.new)
}


marginal_survival <- function(tau, Bh, eo, 
                              surv0, delta0, K,
                              beta, ndx, deg, p, ns=10){
  ## marginal survival estimate
  ntau <- length(tau); nparam <- length(beta)
  n <- length(surv0)
  nbeta0 <- (ndx+deg+1)
  beta0 <- beta[1:nbeta0]
  beta1 <- beta[nbeta0 + 1:p]
  
  Y0 <- seq(1/ns,1,length.out = ns)
  w <- Cost_bbase(Y0, 0, 1, ndx, deg)
  K0 <- seq(11,ntau,10)
  K <- sapply(surv0,function(x) which.max(x<=tau))
  Ytilde <- lapply(1:ns, function(x) c(rep(0,x), 1))
  
  eBhb <- exp(Bh %*% beta0)
  eXb <- exp(X0 %*% beta1)
  rate <- lapply(1:n, function(i) eBhb[1:K[i]] * eo[[i]])
  
  o <- outer(Y0, c(tau,2), FUN=function(x,y) pmin(x,y))
  o1 <- o[,1:(ntau-1)]
  o2 <- o[,3:(ntau+1)]
  o <- apply(log((o2-o1)/2),1,function(x) x[x>-30])
  eo0 <- lapply(1:ns, function(i) c(pmin(tau[2],Y0[i]),exp(o[[i]]))) # exp(o_ik)
  
  rate0 <- lapply(1:ns, function(i) eBhb[1:K0[i]] * eo0[[i]])
  S0 <- exp(-eXb %*% sapply(rate0, sum)) 
  S <- exp(-eXb * sapply(rate, sum)) # marginal survival distribution at 0.1,0.2,...,1
  dS0 <- lapply(1:n,function(i)
    sapply(1:ns, function(j){
      if(S[i]>S0[i,j] & delta0[i]==0){
        round(-S0[i,j] / S[i] *
                c(colSums(matrix(rate0[[j]][K[i]:K0[j]] * 
                                   Bh[K[i]:K0[j],],ncol=nbeta0)),
                  X0[i,]),3)
      }else{rep(0,nparam)}
    }
    )) # deriv of P(s>l_ik|s>Ti,X) wrt beta, conditional on survival time and X
  return(list(surv0 = S0, surv = S,
              deriv_surv0 = dS0))
}


conditional_survival <- function(marg_surv, n, nparam,nbeta0, X0, id0,
                                 surv0, delta0, ns=10){
  ## conditional survival estimate
  
  S0 <- marg_surv$surv0
  S <- marg_surv$surv
  dS0 <- marg_surv$deriv_surv0
  id.cen0 <- which(surv0 < 1 & delta0 == 0)
  surv00 <- surv01 <- seq(1/ns, 1, length.out = ns); surv00[ns] <- 1.1
  prob.cond <- surv0.cond <- matrix(-99, ns+1, n)
  # dprob.cond <- array(0, dim = c(n, nparam, ns+1))
  for(x in id.cen0){
    # calculate derivative Pim(beta)
    idx <- x
    scond <- (surv00 > surv0[x]) *  S0[x,]
    idx0 <- which(scond==0)
    nidx0 <- length(idx0)
    idx0 <- ifelse(nidx0 == 0, 0, idx0[nidx0])
    scond <- append(scond, S[x], after=idx0) / S[x]
    pcond <- c(-diff(scond), scond[ns+1])
    pcond[pcond<0] <- 0
    dcond <- append(surv01, surv0[x], after=idx0)
    dcond <- dcond + c(diff(dcond) / 2, .1)
    prob.cond[,idx] <- pcond
    surv0.cond[,idx] <- dcond
    # calculate derivative d Pim(beta) / d betap
    # dscond <- dS0[[x]] #- dS[[x]]
    # dscond <- unname(data.frame(append(data.frame(dscond), 
    #                                    list(X0=c(rep(0,nbeta0),-X0[x,])), idx0)))
    # dpcond <- cbind(t(-diff(t(dscond))), dscond[,ns+1])
    # dprob.cond[idx,,] <- round(dpcond, 3)
  }
  
  # colnames(surv0.cond) <- colnames(prob.cond) <- paste0('id_',1:n)
  
  # rownames(surv0.cond) = colnames(prob.cond) = id0
  # dprob.cond <- aperm(dprob.cond, c(2,3,1))
  surv0.prob.cond <- list(surv = t(surv0.cond), 
                          prob = t(prob.cond)#, dprob = dprob.cond
                          )
  return(surv0.prob.cond)
}


variance_survival <- function(n, K, Bh, Ytilde, eo, P, beta,
                              surv0, delta0, ndx, deg, p){
  # calculate sandwich variance estimator
  
  beta0 <- beta[1:(ndx+deg+1)]
  beta1<- beta[ndx+deg+1 + 1:p]
  eBhb <- exp(Bh %*% beta0)
  eXb <- exp(X0 %*% beta1)
  
  rate <- lapply(1:n, function(i) eBhb[1:K[i]] * eXb[i] * eo[[i]])
  
  gradn <- round(sapply(1:n, function(i) 
    Reduce('+',lapply(1:K[i],function(j) 
      c(Bh[j,],X0[i,])  * (Ytilde[[i]][j] - rate[[i]][j])))),3) #- c(P %*% beta)
  
  hess <- round(Reduce('+',lapply(1:n, function(i) 
    Reduce('+',lapply(1:K[i],function(j) 
      tcrossprod(c(Bh[j,],X0[i,])) * rate[[i]][j])))) + P,3)
  
  # meat <- Reduce('+',lapply(1:n,function(x) gradn[,x]%*%t(gradn[,x])))
  
  return(list(gradn = gradn, hess = hess))
  
}
