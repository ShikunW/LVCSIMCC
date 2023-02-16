
Cost_bbase <-
  function(x, xl, xr, ndx, deg){
    ## Input:
    ## x   = abcissae of data
    ## xl  = left boundary
    ## xr  = right boundary
    ## ndx = number of internal knots -1
    ##       or number of internal intervals
    ## deg = degree of the splines
    
    ## Output:
    ## B = matrix with the B-spline basis
    ndx <- ndx + 1
    ## distance between knots
    dx <- (xr - xl) / ndx
    ## One needs (ndx+1) internal knots and 
    ## deg knots on both right and left side
    ## in order to joint all the B-splines
    knots <- seq(xl - deg * dx, xr + deg * dx, by = dx)
    #knots <- c(seq(xl - deg * dx, xl + deg * dx, by = dx),seq(xr - deg * dx, xr + deg * dx, by = dx))
    ## Truncated deg-th power functions
    ## equally-spaced at given knots along x
    P <- outer(x, knots, Cost_tpower, deg)
    ## number of B-splines which equal to the number of knots
    n <- dim(P)[2]
    ## in the numerator we have the matrix
    ## of deg+1 differences for each knots
    ## this matrix is rescaled by 
    ## (deg! * dx^deg) == (gamma(deg + 1) * dx ^ deg)
    D <- diff(diag(n), diff = deg + 1) /
      (gamma(deg + 1) * dx ^ deg)
    ## the matrix of differences is used to compute B-splines
    ## as differences of truncated power functions
    ## in P %*% t(D)
    ## the last factor is (-1) ^ (deg + 1)
    B <- (-1) ^ (deg + 1) * P %*% t(D)
    B
  }
Cost_tpower <-
  function(x, t, p){
    ## Input:
    ## x = abcissae of data
    (x - t) ^ p * (x > t)
    ## (x-t)^p gives the curve
    ## (x>t) is an indicator function; it is 1 when x>t
    ## and 0 when x<=t, i.e. before each knot
  }

Cost_pbase <-
  function(x, xl, xr, ndx, deg){
    ## Input:
    ## x   = abcissae of data
    ## xl  = left boundary
    ## xr  = right boundary
    ## ndx = number of internal knots -1
    ##       or number of internal intervals
    ## deg = degree of the splines
    
    ## Output:
    ## B = matrix with the B-spline basis
    P <- c() 
    if(ndx > 0){
      dx <- (xr - xl) / (ndx+1)
      knots <- seq(xl + dx, xr - dx, by = dx)
      P <- outer(x, knots, Cost_tpower, deg)
    }
    P0 <- matrix(sapply(0:deg, function(xx) x^xx),ncol=deg+1)
    B <- cbind(P0, P)
    B
  }


Cost_pbase_deriv <-
  function(x, xl, xr, ndx, deg){
    ## Input:
    ## x   = abcissae of data
    ## xl  = left boundary
    ## xr  = right boundary
    ## ndx = number of internal knots -1
    ##       or number of internal intervals
    ## deg = degree of the splines
    
    ## Output:
    ## B = matrix with the B-spline basis
    P <- c() 
    if(ndx > 0){
      dx <- (xr - xl) / (ndx+1)
      knots <- seq(xl + dx, xr - dx, by = dx)
      P <- outer(x, knots, Cost_tpower, deg-1) * deg
    }
    P0 <- matrix(cbind(rep(0,length(x)),
                       sapply(1:deg, function(xx) xx * x^(xx-1))),ncol=deg+1)
    B <- cbind(P0, P)
    B
  }
Cost_tpower <-
  function(x, t, p){
    ## Input:
    ## x = abcissae of data
    (x - t) ^ p * (x > t)
    ## (x-t)^p gives the curve
    ## (x>t) is an indicator function; it is 1 when x>t
    ## and 0 when x<=t, i.e. before each knot
  }

Convergence <- function( new, old ){
  p <- 2
  n.predictor <- length(new)
  ind <- rep(0, n.predictor)
  for (i in 1:n.predictor) {
    if (new[i] < abs(0.01)){
      ind[i] <- as.numeric( abs(new[i]-old[i]) < 10^(-p) )
    } else{
      ind[i] <- as.numeric( abs(new[i]-old[i]) / (abs(old[i])+10^(-5)) < 10^(-p) )
    } 
  }
  det <- FALSE
  #cat("converge number = ", sum(ind), "\n")
  if(sum(ind) == n.predictor ) det <- TRUE
  det
}
