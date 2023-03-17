evaluate_SLGP_grid <- function(epsilon, Xx, name_index, lengthscale=1,
                               par_basis_functions, n.approxInt=101, 
                               discretization.size=101){
  u <- seq(0, 1,, n.approxInt)
  if(is.vector(Xx)){
    Xx <- matrix(Xx)
    colnames(Xx) <- name_index
  }
  X <- expand.grid(u, seq(nrow(Xx)))
  X <- data.frame(cbind(X[,1], Xx[X[, 2], ]))
  colnames(X) <- c("t", name_index)
  
  storeIntermediate <-storeIntermediate.likelihood(samples=X, 
                                                   par_basis_functions=par_basis_functions, 
                                                   n.approxInt = n.approxInt,
                                                   discretization.size = discretization.size,
                                                   name_index=name_index,
                                                   lengthscale=lengthscale)
  matweight.approxInt<- storeIntermediate$matweight.approxInt
  factormat.normcte<- storeIntermediate$factormat.normcte
  
  factormat.samp<- storeIntermediate$factormat.samp
  map.samp.index<- storeIntermediate$map.samp.index
  map.loc.index<- storeIntermediate$map.loc.index
  df <-data.frame(X)
  colnamesdf <- c(colnames(df))
  if(!is.matrix(epsilon)){
    epsilon_all <- matrix(epsilon, ncol=par_basis_functions$order_tot)
  }else{
    epsilon_all <- epsilon
  }
  for(j in seq(nrow(epsilon_all))){
    epsilon <- epsilon_all[j,]
    df$newfun <-  (factormat.samp %*% epsilon)[map.samp.index, ]
    df <- df%>%
      group_by_at(name_index)%>%
      mutate(newfun = newfun-max(newfun))%>%
      mutate(newfun = exp(newfun))%>%
      mutate(newfuncdf = cumsum(newfun*matweight.approxInt[1:n.approxInt,1]))%>%
      mutate(newfuncdf = newfuncdf-min(newfuncdf))%>%
      mutate(newfun=newfun/max(newfuncdf),
             newfuncdf=newfuncdf/max(newfuncdf))%>%
      ungroup()%>%
      data.frame()
    colnamesdf <-c(colnamesdf, give_names_columns(j, j, "pdf_"), give_names_columns(j, j, "cdf_"))
    colnames(df) <- colnamesdf
  }
  if(nrow(epsilon_all)==1){
    colnames(df) <- c(colnames(X), "pdf", "cdf")
  }
  return(df)
}
