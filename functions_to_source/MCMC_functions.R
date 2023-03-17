library(LaplacesDemon) # For pCN


Model <- function(parm, Data){
  #Calculate log posterior LP
  epsilon <- parm
  storeIntermediate<-Data$storeIntermediate
  starting_lengthscale<-Data$starting_lengthscale
  starting_sigma<-Data$starting_sigma
  likelihood_compute <- - compute.negloglikelihood(epsilon = epsilon,
                                                   factormat.samp = storeIntermediate$factormat.samp, 
                                                   factormat.normcte = storeIntermediate$factormat.normcte,
                                                   matweight.approxInt = storeIntermediate$matweight.approxInt,
                                                   map.loc.index = storeIntermediate$map.loc.index,
                                                   map.samp.index = storeIntermediate$map.samp.index)
  prior_compute <- - compute.neglogpriors(epsilon = epsilon,
                                          sigma=starting_sigma,
                                          lengthscale = starting_lengthscale)
  LP = likelihood_compute + prior_compute
  #Return
  Modelout <- list(LP = LP, Dev = 0, Monitor = LP,
                   yhat = c(0), parm = parm)
  return(Modelout)
}
Model <- cmpfun(Model)

run_pCN <- function(starting_point=NULL, samples, Thinning, Iterations, beta,
                    par_basis_functions, name_index, starting_lengthscale, starting_sigma,
                    n.approxInt = 101, 
                    discretization.size = 101,
                    option=c("pCN", "MALA")){
  storeIntermediate <- storeIntermediate.likelihood(samples=samples, 
                                                    par_basis_functions=par_basis_functions, 
                                                    n.approxInt = n.approxInt,
                                                    discretization.size = discretization.size,
                                                    name_index=name_index,
                                                    lengthscale=starting_lengthscale)
  
  order_tot <- par_basis_functions$order_tot
  if(is.null(starting_point)){
    Initial.Values <- starting_sigma*rnorm(order_tot)
  }else{
    Initial.Values <- starting_point
  }
  parm.names <- paste0("epsilon", seq(order_tot))
  MyData <- list(J = order_tot, N=nrow(samples), X = as.matrix(samples), mon.names = "LP",
                 parm.names = parm.names, storeIntermediate=storeIntermediate,
                 starting_lengthscale=starting_lengthscale,
                 starting_sigma=starting_sigma)
  if(option=="pCN"){
    Fit <- LaplacesDemon(Model, Data = MyData, Initial.Values,
                         Covar = NULL, Iterations = Iterations, Status = 1e3,
                         Thinning = Thinning,
                         Algorithm = "pCN", Specs = list(beta = beta))
    
    # Initial.Values <- as.initial.values(Fit)
    # Fit <- LaplacesDemon(Model, Data=MyData, Initial.Values,
    #                      Covar=Fit$Covar, Iterations=1e+05, Status=10204, Thinning=1000,
    #                      Algorithm="RWM")
  }
  else if(option=="MALA"){
    Fit <- LaplacesDemon(Model, Data = MyData, Initial.Values,
                         Covar = NULL, Iterations = Iterations, Status = 1e1,
                         Thinning = Thinning,
                         Algorithm="MALA", Specs=list(A=1e7, alpha.star=0.574, gamma=1,
                                                      delta=1, epsilon=c(1e-6,1e-7)))
  }else{
    Fit <- LaplacesDemon(Model, Data = MyData, Initial.Values,
                         Covar = NULL, , Iterations=Iterations, 
                         Status=1000, Thinning=Thinning,
                         Algorithm="AM", Specs=list(Adaptive=500, Periodicity=10))
  }
  # Consort(Fit)
  
  return(Fit)
}