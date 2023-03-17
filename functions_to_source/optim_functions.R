
compute_MAP_without_hyperpar <- function(samples, name_index,
                                         par_basis_functions,
                                         starting_point,
                                         starting_sigma =1,
                                         starting_lengthscale = 1,
                                         print_level=0){
  order_tot <- par_basis_functions$order_tot
  dimension <- par_basis_functions$dimension
  starting_lengthscale <- c(starting_lengthscale)
  
  if(length(starting_lengthscale)==1){
    starting_lengthscale <- rep(starting_lengthscale, dimension)
  }
  
  
  ordered_columns <- colnames(samples)
  ordered_columns <- c(ordered_columns[!(ordered_columns %in% name_index)], name_index)
  samples <- samples[, ordered_columns]
  # samples_scaled <- data.frame(t(t(samples)/starting_lengthscale))
  
  storeIntermediate <-storeIntermediate.likelihood(samples, 
                                                   par_basis_functions, 
                                                   n.approxInt = 101,
                                                   discretization.size = 101,
                                                   name_index=name_index,
                                                   lengthscale=starting_lengthscale)
  funoptim <- function(x){
    epsilon <- x
    value <- compute.neglogpriors(epsilon = epsilon,
                                  sigma=starting_sigma,
                                  lengthscale = starting_lengthscale) +
      compute.negloglikelihood(epsilon = epsilon,
                               factormat.samp = storeIntermediate$factormat.samp, 
                               factormat.normcte = storeIntermediate$factormat.normcte,
                               matweight.approxInt = storeIntermediate$matweight.approxInt,
                               map.loc.index = storeIntermediate$map.loc.index,
                               map.samp.index = storeIntermediate$map.samp.index)
    return(value)
  }
  gradfunoptim <- function(x){
    epsilon <- x
    value_prior <- numDeriv::grad(function(y){ 
      compute.neglogpriors(epsilon = y,
                           sigma=starting_sigma,
                           lengthscale = starting_lengthscale)}, x=epsilon)
    value <- value_prior +
      compute.grad.negloglikelihood(epsilon = epsilon,
                                    factormat.samp = storeIntermediate$factormat.samp, 
                                    factormat.normcte = storeIntermediate$factormat.normcte,
                                    matweight.approxInt = storeIntermediate$matweight.approxInt,
                                    map.loc.index = storeIntermediate$map.loc.index,
                                    map.samp.index = storeIntermediate$map.samp.index)
    return(value)
  }
  if(FALSE){
    #Gradient seems ok!
    library(numDeriv)
    for(i in seq(3)){
      x <- rnorm(order_tot, mean=0, sd=starting_sigma)
      print(range(abs(gradfunoptim(x)-numDeriv::grad(funoptim, x))))
    }
  }
  require(nloptr)
  res_optim <- nloptr::nloptr(starting_point, eval_f=funoptim,
                      eval_grad_f = gradfunoptim,
                      lb=c(-100+0*starting_point),
                      ub=c(100+0*starting_point),
                      opts=list(algorithm="NLOPT_LD_LBFGS",
                                xtol_rel = 1e-04, ftol_rel = 1e-04,
                                maxeval=10000,
                                print_level=0))
  # res_optim <- nloptr(starting_point, eval_f=funoptim,
  #                     lb=c(-Inf+0*starting_point),
  #                     ub=c(Inf+0*starting_point),
  #                     opts=list(algorithm="NLOPT_LN_COBYLA",
  #                               xtol_rel = 1e-04, ftol_rel = 1e-04,
  #                               maxeval=10000,
  # #                               print_level=0))
  # res_optim <- optim(par=starting_point, fn=funoptim,
  #                    gr=gradfunoptim, method = "L-BFSG",
  #                    control=list(trace=print_level))
  epsilon_opt <- res_optim$solution
  value<-res_optim$objective
  # value<-res_optim$value
  # epsilon_opt <- res_optim$par
  sigma_opt <-starting_sigma
  lengthscale_opt <- starting_lengthscale
  
  return(list(epsilon=epsilon_opt,
              sigma=sigma_opt,
              lengthscale=lengthscale_opt,
              value=value))
}

