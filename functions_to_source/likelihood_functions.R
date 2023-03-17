require(compiler)
storeIntermediate.likelihood <- function(samples, par_basis_functions, lengthscale=1,
                                         n.approxInt = 101,
                                         discretization.size = 101,
                                         name_index){
  require(dplyr)
  #Pre-treat samples
  if(!is.null(discretization.size)){
    samples <- round(samples*(discretization.size-1))/(discretization.size-1)
  }
  loc <- as.matrix(dplyr::select(samples, all_of(name_index)))
  #Map to a grid to reduce memory required : 
  unique.samp <- unique(samples)
  map.samp.index <- match(data.frame(t(samples)), data.frame(t(unique.samp)))
  
  unique.loc <- unique(loc)
  unique.loc.index <- seq(nrow(unique.loc))
  map.loc.index <- match(data.frame(t(loc)), data.frame(t(unique.loc)))
  
  #Approximate integrals
  n.distinct.index <- length(unique.loc.index)
  n.samples <- nrow(samples)
  dim.index <- ncol(loc)
  
  #To approximate the normalization constant
  weight.approxInt <- c(1, rep(c(4, 2), (n.approxInt-3)/2), 4,1) / (3*(n.approxInt-1))
  matweight.approxInt <- matrix(0, nrow=n.approxInt*n.distinct.index, ncol = n.distinct.index)
  for(i in seq(n.distinct.index)){
    matweight.approxInt[(i-1)*n.approxInt+1:n.approxInt, i] <- weight.approxInt
  }
  u.approxInt <- seq(0, 1,, n.approxInt)
  vect.approxInt <- expand.grid(u.approxInt, unique.loc.index)
  vect.approxInt <- cbind(vect.approxInt[, 1], unique.loc[vect.approxInt[, 2], ])
  colnames(vect.approxInt) <- c("t", name_index)
  
  factormat.normcte <- evaluate_basis_functions(X=vect.approxInt, 
                                                par_basis_functions=par_basis_functions,
                                                lengthscale=lengthscale)
  factormat.samp  <- evaluate_basis_functions(X=unique.samp, 
                                              par_basis_functions=par_basis_functions,
                                              lengthscale=lengthscale)
  
  return(list(factormat.samp = factormat.samp, 
              factormat.normcte = factormat.normcte, 
              matweight.approxInt = matweight.approxInt,
              map.loc.index=map.loc.index,
              map.samp.index=map.samp.index))
}
storeIntermediate.likelihood <- cmpfun(storeIntermediate.likelihood)

compute.likelihood <- function(epsilon, 
                               factormat.samp, 
                               factormat.normcte,
                               matweight.approxInt,
                               map.loc.index,
                               map.samp.index){
  #For numerical stability
  funvalInt <- c(factormat.normcte %*% epsilon)
  funvalSamp <- c(factormat.samp %*% epsilon)
  funvalSamp <- funvalSamp[map.samp.index] 
  nx <- ncol(matweight.approxInt)
  n.approxInt <- nrow(matweight.approxInt)/nx
  group_funvalInt <- sapply(seq(nx), function(x){rep(x, n.approxInt)})
  stability_value <- unname(tapply(funvalInt, group_funvalInt, max))
  funvalInt <- funvalInt-stability_value[group_funvalInt]
  funvalSamp <- funvalSamp -stability_value[map.loc.index]
  
  norm.cte <- c(exp(funvalInt) %*% matweight.approxInt)
  likelihood <- c(exp(funvalSamp) / norm.cte[map.loc.index])
  return(likelihood)
}
compute.likelihood <- cmpfun(compute.likelihood)

compute.negloglikelihood <- function(epsilon, 
                                     factormat.samp, 
                                     factormat.normcte,
                                     matweight.approxInt,
                                     map.loc.index,
                                     map.samp.index){
  #For numerical stability
  funvalInt <- c(factormat.normcte %*% epsilon)
  funvalSamp <- c(factormat.samp %*% epsilon)
  funvalSamp <- funvalSamp[map.samp.index] 
  nx <- ncol(matweight.approxInt)
  n.approxInt <- nrow(matweight.approxInt)/nx
  group_funvalInt <- sapply(seq(nx), function(x){rep(x, n.approxInt)})
  stability_value <- unname(tapply(funvalInt, group_funvalInt, max))
  funvalInt <- funvalInt-stability_value[group_funvalInt]
  funvalSamp <- funvalSamp -stability_value[map.loc.index]
  
  norm.cte <- c(exp(funvalInt) %*% matweight.approxInt)
  negloglikelihood <- sum(-funvalSamp + log(norm.cte[map.loc.index]))
  return(negloglikelihood)
}
compute.negloglikelihood <- cmpfun(compute.negloglikelihood)


compute.grad.negloglikelihood <- function(epsilon, 
                                          factormat.samp, 
                                          factormat.normcte,
                                          matweight.approxInt,
                                          map.loc.index,
                                          map.samp.index){
  #For numerical stability
  funvalInt <- c(factormat.normcte %*% epsilon)
  funvalSamp <- c(factormat.samp %*% epsilon)
  funvalSamp <- funvalSamp[map.samp.index] 
  nx <- ncol(matweight.approxInt)
  n.approxInt <- nrow(matweight.approxInt)/nx
  group_funvalInt <- sapply(seq(nx), function(x){rep(x, n.approxInt)})
  stability_value <- unname(tapply(funvalInt, group_funvalInt, max))
  funvalInt <- funvalInt-stability_value[group_funvalInt]
  funvalSamp <- funvalSamp -stability_value[map.loc.index]
  
  # First term of the gradient
  grad1 <- - colSums(factormat.samp[map.samp.index, ])
  
  # Second term
  norm.cte <- c(exp(funvalInt) %*% matweight.approxInt)
  term.density.grad <- c(exp(funvalInt) / norm.cte[group_funvalInt])
  temp <- factormat.normcte*term.density.grad
  grad2 <- colSums(apply(temp, 2, function(y){
    tapply(y, group_funvalInt, function(x){
      sum(x*matweight.approxInt[1:n.approxInt, 1])
    })
  })[map.loc.index, ])
  
  grad <- grad1+grad2
  return(grad)
}
compute.grad.negloglikelihood <- cmpfun(compute.grad.negloglikelihood)


# numericalgradient <- grad(function(x){
#   return(-sum(log(compute.likelihood(epsilon=x,
#                                      factormat.samp=factormat.samp,
#                                      factormat.normcte=factormat.normcte,
#                                      matweight.approxInt=matweight.approxInt,
#                                      map.loc.index=map.loc.index,
#                                      map.samp.index=map.samp.index))))
# }, x=epsilon)
# analyticalgradient <- compute.grad.likelihood(epsilon=epsilon,
#                                  factormat.samp=factormat.samp,
#                                  factormat.normcte=factormat.normcte,
#                                  matweight.approxInt=matweight.approxInt,
#                                  map.loc.index=map.loc.index,
#                                  map.samp.index=map.samp.index)
# range(numericalgradient - analyticalgradient)
