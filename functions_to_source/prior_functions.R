neglogpriorepsilon <- function(epsilon, sigma){
  #For known sigma
  fun <- -sum(dnorm(epsilon, mean=0, sd=sigma, log=TRUE))
  return(list(fun = fun))
}

neglogpriorsigma <- function(sigma,
                             differentiate=FALSE){
  
  #To adapt, here I decided for a Gamma distributions
  shape <- 3
  scale <- 1/2
  fun <- sigma/scale - (shape-1)*log(sigma)+shape*log(scale)+lgamma(shape)
  fun <- sum(fun)
  if(differentiate){
    grad <- 1/scale - (shape-1)/sigma
    hess <- diag((shape-1)/sigma^2)
    return(list(fun = fun, grad = grad, hess = hess))
  }
  return(list(fun = fun))
}

neglogpriorlengthscale <- function(lengthscale, differentiate=FALSE, verbose=0){
  
  require(invgamma)
  shape <- 8
  mode_at <- 0.1
  rate <- mode_at*(shape+1)
  if(verbose==1){
    # u <- seq(0, 1, 0.001)
    # y <- dinvgamma(u, shape=shape, rate=rate)
    meanval <- rate/(shape-1)
    modeval <- rate/(shape+1)
    q1 <- invgamma::qinvgamma(shape=shape, rate=rate, 0.05/2)
    q2 <- invgamma::qinvgamma(shape=shape, rate=rate, 1-0.05/2)
    # plot(u, y, "l")
    # abline(v=c(q1, q2), lty=2)
    # abline(v=c(meanval, modeval), lty=3, col="grey")
    print(paste0("With shape: ", shape, " and rate: ", rate, 
                 ", the mode will be: ", round(modeval, 2),
                 ", the mean will be: ", round(meanval, 2),
                 " and 95% of the values will be in [",
                 round(q1, 2),
                 ", ",
                 round(q2, 2),
                 "]."
    ))
  }
  func <- function(lengthscale){
    return(-sum(invgamma::dinvgamma(lengthscale, shape=shape, rate=rate, log=TRUE)))
  }
  fun <- func(lengthscale)
  if(differentiate){
    require(numDeriv)
    grad <- grad(func, x=lengthscale)
    hess <- hessian(func, x=lengthscale)
    return(list(fun = fun, grad = grad, hess = hess))
  }
  return(list(fun = fun))
}


compute.neglogpriors <- function(epsilon, 
                                 sigma,
                                 lengthscale){
  neglogprior <- neglogpriorepsilon(epsilon, sigma)$fun+
    neglogpriorsigma(sigma)$fun+
    neglogpriorlengthscale(lengthscale)$fun
  return(neglogprior)
}
compute.neglogpriors<- cmpfun(compute.neglogpriors)
