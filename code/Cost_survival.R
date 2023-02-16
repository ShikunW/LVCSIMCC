library(survival)
library(MASS)
library(parallel)
# library(glmnet)
# dat=Cost_generate_data(n = 500, prop.zero = 0,
#                          seed = 123, test = F, fig = T, dist = 1)
# dat.surv <- dat[!duplicated(dat$id),]
# surv0 <- dat.surv$surv
# delta0 <- dat.surv$delta
# X0 <- as.matrix(dat.surv[,substr(names(dat.surv),0,1)=='x'])

Cost_survival0 <- function(surv0, delta0, X0){
  # Initial value for theta_3
  # deg <- 1; ndx <- 1
  n <- length(surv0)
  ns <- 10
  colnames(X0) <- paste0('x',1:ncol(X0))
  dat.surv <- data.frame(id = 1:n, surv0 = surv0, delta0 = delta0, X0 = X0)
  survformula <- formula(paste0("Surv(surv0,delta0) ~ X0.x1",paste(paste0('+X0.x',2:ncol(X0)), collapse = '')))
  model <- coxph(survformula, data = dat.surv)
  
  s0 <- seq(.1,1,length.out = ns)
  CondProb <- matrix(-99,n,ns); CondSurv <- matrix(1.01,n,ns)
  for(i in which(surv0 < 1 & delta0 == 0)){
    idx.this <- which(s0>surv0[i])
    s0.this <- c(surv0[i], s0[idx.this])
    surv.this <- s0.this[1:(length(s0.this))] 
    CondSurv[i,idx.this] <- s0.this[1:(length(s0.this)-1)] + diff(s0.this)/2
    newdata <- data.frame(surv0 = surv.this,
                          delta0 = rep(1, length(surv.this)))
    newdata <- cbind(newdata, matrix(X0[rep(i, each=length(surv.this)),],ncol = ncol(X0)))
    names(newdata)[3:ncol(newdata)] <- paste0('X0.x',names(newdata)[3:ncol(newdata)])
    
    pred.this <- predict(model,newdata=newdata, type='survival')
    pred.this <- pred.this / pred.this[1]
    CondProb.this <- -diff(pred.this)
    CondProb[i,idx.this] <- CondProb.this
  }
  
  CondProb <- cbind(CondProb,1-apply(CondProb,1,sum))
  CondSurv <- cbind(CondSurv, rep(1.1,n))
  surv.out <- list(surv = CondSurv,  prob = CondProb, coefficient = coef(model))
  return(surv.out)
}

Cost_survival <- function(surv0, delta0, X0,
                          ndx = 10, deg = 3, lambda = .001){
  # Initial value for theta_3
  # deg <- 1; ndx <- 1
  theta3.new <- matrix(coef(coxph(Surv(surv0,delta0) ~ X0)))
  theta3.new <- theta3.new / sqrt(sum(theta3.new^2))
  ntheta3 <- ndx + deg + 1
  beta.new <- rep(0,2*ntheta3 - 1)
  
  # Poissonization
  KK <- 101
  tau0 <- c(seq(0,1,length.out = KK), 1.01)
  Bp0 <- Cost_bbase(tau0, 0, 1, ndx, deg)
  n <- length(surv0)
  sl <- 0
  sr <- 1
  smax <- sr + 0.01 * (sr - sl)
  smin <- sl - 0.01 * (sr - sl)
  # penalty matrix
  DtD <- diag(c(rep(0,deg),rep(1,ndx)))
  p <- ncol(X0)
  nd <- 103 # M+2
  K <- colSums(matrix(kronecker(surv0, tau0, '>'), nrow = KK+1)) + 1
  o <- matrix(0,nrow = n, ncol = KK)
  Y <- lapply(1:n, function(i) c(rep(0,K[i]-1),delta0[i]))
  yy <- Reduce('c', Y)
  o[,1] <- rep(tau0[2], n)
  o[,2:KK] <- (pmin(matrix(tau0[3:(KK+1)], nrow=n, ncol=KK-1, byrow=TRUE),
                    t(matrix(surv0, nrow=KK-1, ncol=n, byrow=TRUE))) -
                 pmin(matrix(tau0[1:(KK-1)], nrow=n, ncol=KK-1, byrow=TRUE),
                      t(matrix(surv0, nrow=KK-1, ncol=n, byrow=TRUE)))) / 2
  lo <- log(o)
  lo[lo==-Inf] <- 0
  oo <- Reduce('c', sapply(1:n, function(i) lo[i,1:K[i]]))
  ns <- 10
  s0 <- seq(.1,1,length.out = ns)
  o1 <- matrix(0,nrow = ns, ncol = KK)
  o1[,1] <- rep(tau0[2], ns)
  o1[,2:KK] <- (pmin(matrix(tau0[3:(KK+1)], nrow=ns, ncol=KK-1, byrow=TRUE),
                     t(matrix(s0, nrow=KK-1, ncol=ns, byrow=TRUE))) -
                  pmin(matrix(tau0[1:(KK-1)], nrow=ns, ncol=KK-1, byrow=TRUE),
                       t(matrix(s0, nrow=KK-1, ncol=ns, byrow=TRUE)))) / 2
  
  iterate <- 0; MAX.IT <- 20
  repeat{
    theta3 <- theta3.new
    beta <- beta.new
    Xtheta3 <- c(X0 %*% theta3)
    B3 <- Cost_bbase(Xtheta3, 0, 1, ndx, deg)[,-1]
    BB <- lapply(1:n, function(i) cbind(B3[rep(i,K[i]),],Bp0[1:K[i],]))
    Z <- Reduce('rbind',BB)
    
    p.fac <- rep(c(rep(0,deg),rep(1,ndx)), 2)
    
    if( ndx > 0 ){
      model.1 <- glmnet(y = yy, x = Z, offset = oo, intercept = F,
                        family = poisson(link = "log"), alpha=0, lambda=0,
                        penalty.factor = p.fac)
      beta.new <- as.array(coef(model.1)[-1])
    }else{
    model.1 <- glm(yy ~ -1+I(Z), offset = oo,
                   family = poisson(link = "log"))

    beta.new <- as.array(coef(model.1))
    }
    
    beta.new[is.na(beta.new)] <- 0
    beta1.new <- beta.new[1:(ntheta3-1)]
    beta2.new <- beta.new[1:(ntheta3) + ntheta3-1]
    Bp1b <- Cost_bbase_deriv(Xtheta3, 0, 1, ndx, deg)[,-1] %*% beta1.new
    Jt1 <- rbind(-theta3[2:p] / sqrt(1-sum(theta3[2:p]^2)), diag(rep(1,p-1)))
    ZZ <- c(t(Jt1) %*% t(X0)) * Bp1b
    Z1 <- Reduce('c', sapply(1:n, function(i) rep(ZZ[i],K[i])))
    Z1 <- matrix(Z1)
    tmp <- Bp1b %*% theta3[2:p]
    tmp <- Reduce('c',sapply(1:n, function(i) tmp[rep(i,K[i])]))
    o2 <- Z %*% beta.new - tmp
    
    model.2 <- glm(yy ~ -1+I(Z1), offset = oo+o2,
                   family = poisson(link = "log"))
    
    theta3.new <- c(coef(model.2))
    theta3.new <- c(1-sum(theta3.new^2),theta3.new)
    
    iterate <- iterate + 1
    print(iterate)
    if( Convergence(c(theta3.new, beta.new), c(theta3, beta)) | iterate > MAX.IT){
      converge <- 1
      break
    }
  }
  
  Xtheta3 <- c(X0 %*% theta3.new)
  B3 <- Cost_bbase(Xtheta3, 0, 1, ndx, deg)[,-1]
  BB <- lapply(1:n, function(i) cbind(B3[rep(i,nd-1),],Bp0))
  
  W <- Cost_bbase(surv0, xl=0, xr=1, ndx, deg)
  
  X1 <- matrix(c(0,0,.5,.5,0,1,0,1),ncol=2)
  X1theta3 <- c(X1 %*% theta3.new)
  B13 <- Cost_bbase(X1theta3, 0, 1, ndx, deg)[,-1]
  BB1 <- lapply(1:nrow(X1), function(i) cbind(B13[rep(i,nd-1),],Bp0))
  
  exp(-as.vector(exp(BB1[[1]] %*% beta.new) / 100)[1:KK])
  
  ts0.idx <- which(tau0 %in% s0)
  CondProb <- matrix(0,n,ns); CondSurv <- matrix(1.01,n,ns)
  for(i in 1:n){
    CumHaz0 <- as.vector(exp(BB[[i]] %*% beta.new))[1:KK]
    CumHaz1 <- sum(o[i,] * CumHaz0)
    idx.this <- s0>surv0[i]
    s0.this <- c(surv0[i], s0[idx.this])
    CondSurv[i,idx.this] <- s0.this[1:(length(s0.this)-1)] + diff(s0.this)/2
    CumHaz2 <- colSums(matrix(t(o1[idx.this,]),nrow=KK) * CumHaz0)
    CondProb.this <- -diff(exp(-c(0, CumHaz2-CumHaz1)))
    CondProb[i,idx.this] <- CondProb.this
  }
  CondProb <- cbind(CondProb,1-apply(CondProb,1,sum))
  CondSurv <- cbind(CondSurv, rep(1.1,n))
  surv.out <- list(surv = CondSurv,  prob = CondProb)
  return(surv.out)
}
