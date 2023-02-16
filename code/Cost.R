# library(svcm)
# dat=Cost_generate_data(n = 500, prop.zero = 0,
#                          seed = 123, test = F, fig = T, dist = 1)
# 
# surv <- dat$surv
# delta <- dat$delta
# time <- dat$time
# surv <- dat$surv
# Y <- dat$Y
# id <- dat$id
# X <- as.matrix(dat[,substr(names(dat),0,1)=='x'])
# ndx=4;deg=2
# censor=1
Cost <- function(Y, X, time, surv, id=NULL, delta=NULL, 
                 ndx=0, deg=1, ndx0=0, deg0=1, ndx1=0, deg1=1, 
                 correlation='compound', family='normal',constrain=F){
    #time = dat$time; surv = dat$surv; Y = dat$Y;id = dat$id; X=cbind(dat$x1,dat$x2); delta = dat$delta; censor <- T; ndx=5; deg=2; lambdas=c(0,0,0); coefstart=NULL; theta.init = NULL; control=list(); correlation = 'compound'
    
    if(is.null(id)){
      id <- 1:length(Y)
    }
    if(is.null(delta)){
      delta <- rep(1, length(Y))
    }
    X <- as.matrix(X)
    
    time <- time / max(surv)
    surv <- surv / max(surv)
    Y[is.nan(Y)] <- 0; X[is.nan(X)] <- 0; time[is.nan(time)] <- 0; surv[is.nan(surv)] <- 0;
    id0=id[!duplicated(id)]
    surv0 = surv[!duplicated(id)]
    delta0=delta[!duplicated(id)]
    X0=X[!duplicated(id),]
    surv.out <- Cost_survival0(surv0, delta0, X0, ndx0, deg0)#, ndx0, deg0
    # surv.out <- Cost_survival(surv0, delta0, X0, 
    #                           ndx = 1, deg = 1, lambda = .01)
    
    design <- Cost_design(Y, X, time, surv, id, delta, ndx0, deg0, ndx, deg, 
                          censor, surv.out, family = family)
    
    Y=design$Y; X=design$X; time=design$time; surv=design$surv
    id=design$id; delta=design$delta; B=design$B; B0=design$B0; 
    # dB=design$dB; dB0=design$dB0
    rm(design)

    lambda = 1e3
    correlation <- 'compound'
    CondProb <- surv.out$cond_surv$prob
    MAX.IT=2
    param0 <- Cost_estimate(Y, X, B0, B, id, delta, surv, CondProb, 
                            ndx, deg, ndx0, deg0, ndx1, deg1, 
                            param0, lambda, correlation, varfunc='constant', 
                            censor=F, family = family, constrain=F)
    lambda = 1e2
    MAX.IT=2
    # source('Y:/Cost_project/CostPaper3/08-10-2021/Cost_estimate.R')
    param1 <- Cost_estimate(Y, X, B0, B, id, delta, surv, CondProb, 
                            ndx, deg, ndx0, deg0, ndx1, deg1, 
                            param1, lambda, correlation, varfunc='constant', 
                            censor=T, family = family, constrain=F)
    # lambda = 1e3
    # param2 <- Cost_estimate(Y, X, B0, B, id, delta, surv, CondProb, 
    #                         ndx, deg, ndx0, deg0, ndx1, deg1, 
    #                         param1, lambda, correlation, varfunc='spline', 
    #                         censor=F, family = family, constrain=F)
    # param3 <- Cost_estimate(Y, X, B0, B, id, delta, surv, CondProb, 
    #                         ndx, deg, ndx0, deg0, ndx1, deg1, 
    #                         param2, lambda, correlation, varfunc='spline', 
    #                         censor=T, family = family, constrain=F)
    
    variance0 <- Cost_variance(Y, X, B0, B, id, delta, surv, surv.out, 
                               ndx, deg, ndx0, deg0, ndx1, deg1, param0,
                               lambda, MON, correlation, 
                               varfunc='constant', censor=F, family = family)
    surv.out$var_surv$gradn=rbind(surv.out$var_surv$gradn,matrix(0,2,12010))
    surv.out$var_surv$hess=as.matrix(bdiag(surv.out$var_surv$hess,diag(2)))
    variance1 <- Cost_variance(Y, X, B0, B, id, delta, surv, surv.out, 
                               ndx, deg, ndx0, deg0, ndx1, deg1, param1,
                               lambda, MON, correlation, 
                               varfunc='constant', censor=T, family = family)
    variance2 <- Cost_variance(Y, X, B0, B, id, delta, surv, surv.out, 
                               ndx, deg, ndx0, deg0, ndx1, deg1, param2,
                               lambda, MON, correlation, 
                               varfunc='spline', censor=F, family = family)
    variance3 <- Cost_variance(Y, X, B0, B, id, delta, surv, surv.out, 
                               ndx, deg, ndx0, deg0, ndx1, deg1, param3,
                               lambda, MON, correlation, 
                               varfunc='spline', censor=T, family = family)
    # param1 = variance1=0
    # # stopCluster(cl)
    # ## output
    object <- list(
      param0 = param0, variance0 = variance0,
      param1 = param1, variance1 = variance1,
      # param2 = param2, variance2 = variance2,
      # param3 = param3, variance3 = variance3,
      surv = surv.out$coefficient
    )
    class(object) <- "Cost"
    object
  }


