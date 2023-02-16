library(splines)
library(Matrix)
library(survival)
library(MASS)
library(simstudy)
library(evd)
Cost_generate_data <- function(n = 1000, prop.zero = 0, seed = 123, theta1, theta2,
                               test = F, fig = F, family = 'normal'){
  #n = 500; rho = 0.2; prop.zero = 0; seed = 123; test = F; dist = 1; fig = T; family = 'normal'
  set.seed(seed)
  # generate survival data
  if(test){
    n <- 11*4
    x1 <- rep(c(0,0.5), each=n/2)
    x2 <- rep(c(0,1,0,1), each=n/4)
    death <- cbind(seq(0.1,1.1,length.out = n/4),rep(999,n/4))
    survival <- death <- rbind(death,death,death,death)
  }else{
    x1 <- runif(n)
    x2 <- rbinom(n, 1, 0.5)
    U <- runif(n)
    scale =1/3; shape = 1;
    death <- -(log(U) / (scale * exp(2 * x1 + 1 * x2)))^(1/shape)
    death=ceiling(death*100)/100
    censor <- runif(n,0,2)
    death <- cbind(death, censor) 
    survival <- pmax(round(death * 100) / 100, .01)
  }
  
  
  dat.surv <- apply(survival,1,function(x) min(c(x, 1)))
  dat.delta <- apply(survival,1,function(x)1*(x[1] <= x[2] & x[1] <= 1))
  nn <- as.integer(floor(dat.surv * 100))
  dat.tmp <- data.frame(id = rep(1:n,nn),
                        x1 = rep(x1,nn),
                        x2 = rep(x2,nn),
                        surv = rep(dat.surv,nn),
                        death = rep(death[,1],nn),
                        time = unlist(sapply(nn,function(x) 1:x)) / 100,
                        delta = rep(dat.delta,nn))

  X <- as.matrix(dat.tmp[,c('x1','x2')])
  time <- dat.tmp$time
  death <- dat.tmp$death
  if(family == 'normal'){
    dat.tmp$Y <- (cos(2*pi*time/death) + 1) + ((X %*% theta1)^2)
    Ylts <- (cos(2*pi*time) * (time <= 0.5) - (time > 0.5) + 1) + ((X %*% theta2)^2)
    # dat.tmp$Y <- exp(time/death) + ((X %*% theta1)^1)
    # Ylts <- exp(time) + ((X %*% theta2)^1)
  }else if(family %in% c('gamma','poisson')){
    dat.tmp$Y <- (cos(2*pi*time/death)+2)  + ((X %*% theta1)^2)
    Ylts <- (cos(2*pi*time) * (time <= 0.5)- (time > 0.5)+2) + ((X %*% theta2)^2)
    # dat.tmp$Y <- exp(time/death) + ((X %*% theta1)^2)
    # Ylts <- exp(time) + ((X %*% theta2)^1)
  }else if(family == 'binomial'){
    dat.tmp$Y <- exp(time/death) - 2 + ((X %*% theta1)^1 * time * death)
    Ylts <- exp(time) - 2 + ((X %*% theta2)^1 * time)
  }
  dat.tmp$Y[dat.tmp$death>1] <- Ylts[dat.tmp$death>1]
  if(family == 'poisson'){
	  dat.tmp$Y <- exp(dat.tmp$Y)
	  }else if(family == 'binomial'){
	    dat.tmp$Y <- exp(dat.tmp$Y) / (exp(dat.tmp$Y) + 1)
    }
  if(!test){
    for(i in 1:n){
      idxi <- which(dat.tmp$id==i)
      ni <- length(idxi)
      timei <-dat.tmp$time[idxi]
      Yi <- dat.tmp$Y[idxi]
      {
        if(family == 'normal'){ # normal distribution
          rho=0.2
          Ri <- matrix(rho,ni,ni)
          diag(Ri) <- 1 
          varYi.sqrt <- diag(ni)
          diag(varYi.sqrt) <- .16
          if(family=='poisson') diag(varYi.sqrt) <- 1
          Sigma <- (varYi.sqrt) %*% Ri %*% (varYi.sqrt)
          err <- mvrnorm(1,mu = rep(0,ni),Sigma = Sigma) * (dat.tmp$Y[idxi] > 0)
          dat.tmp$Y[idxi] <- Yi + err
        }else if (family == 'gamma'){ # gamma distribution
          if(dat.tmp$death[idxi[1]] <= 1){
            dat.tmp$Y[idxi] <- genCorGen(1,nvars=ni,params1=Yi,params2=1,
                                         dist='gamma',rho=0.2,corstr='cs')$X # mean=params1, var=param1^2*param2
          }else{
            dat.tmp$Y[idxi] <- genCorGen(1,nvars=ni,params1=Yi,params2=0.8,
                                         dist='gamma',rho=0.1,corstr='cs')$X # mean=params1, var=param1^2*param2
          }
        }else if (family == 'poisson'){ # poisson distribution
          if(dat.tmp$death[idxi[1]] <= 1){
            dat.tmp$Y[idxi] <- genCorGen(1,nvars=ni,params1=Yi,
                                         dist='poisson',rho=0.2,corstr='cs')$X # mean=params1, var=param1^2*param2
          }else{
            dat.tmp$Y[idxi] <- genCorGen(1,nvars=ni,params1=Yi,
                                         dist='poisson',rho=0.1,corstr='cs')$X # mean=params1, var=param1^2*param2
          }
        }else if (family == 'binomial'){ # poisson distribution
          if(dat.tmp$death[idxi[1]] <= 1){
            dat.tmp$Y[idxi] <- genCorGen(1,nvars=ni,params1=Yi,
                                         dist='binary',rho=0.2,corstr='cs')$X # mean=params1, var=param1^2*param2
          }else{
            dat.tmp$Y[idxi] <- genCorGen(1,nvars=ni,params1=Yi,
                                         dist='binary',rho=0.1,corstr='cs')$X # mean=params1, var=param1^2*param2
          }
        }
      }
    }
  }
  dat.tmp$Y <- dat.tmp$Y*rbinom(nrow(dat.tmp), p=1-prop.zero, 1) / (1-prop.zero)
  
  # plot data
  if(fig){
    if(test){
      plot(NULL, xlim = c(0,1), ylim = c(0,ifelse(family=='binomial',1,2)), xlab = 'Time', ylab = 'Y')
      for(i in seq(0.1,1.1,.2)){
        dati <- dat.tmp[dat.tmp$death==i & dat.tmp$x1==0 & dat.tmp$x2==0,]
        lines(dati$time,dati$Y,col = floor(dati$death[1]*10))
        dati <- dat.tmp[dat.tmp$death==i & dat.tmp$x1==1 & dat.tmp$x2==0,]
        lines(dati$time,dati$Y,col = floor(dati$death[1]*10),lty=2)
        dati <- dat.tmp[dat.tmp$death==i & dat.tmp$x1==0 & dat.tmp$x2==1,]
        lines(dati$time,dati$Y,col = floor(dati$death[1]*10),lty=3)
        dati <- dat.tmp[dat.tmp$death==i & dat.tmp$x1==1 & dat.tmp$x2==1,]
        lines(dati$time,dati$Y,col = floor(dati$death[1]*10),lty=4)
      }
    }else{
      plot(NULL, xlim = c(0,1), ylim = c(0,2), xlab = 'Time', ylab = 'Y')
      for(i in 1:100){
        dati <- dat.tmp[dat.tmp$id==i,]
        lines(dati$time,dati$Y,col = floor(dati$death[1]*10))
      }
    }
  }
  return(dat.tmp)
  
}
# dat=Cost_generate_data(n = 500, prop.zero = 0,seed = 123, 
#                        theta1=theta1, theta2=theta2,test = T, fig = T,family='normal')
