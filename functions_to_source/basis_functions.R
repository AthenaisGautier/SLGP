require(kergp)

evaluate_basis_function_FF <- function(freq, offset, coef, X, checknames=FALSE){
  if(is.data.frame(X)){
    X <- as.matrix(X)
  }
  Xfreq <- X%*%t(freq)
  basis_fun <- cos(t(Xfreq) + offset)
  basis_fun <- t(coef*basis_fun)
  return(basis_fun)
}

if(FALSE){
  set.seed(3)
  order_tot <- 100
  init <- initialize_basis_function_FF(type="RFF", order_tot=order_tot, dimension=2, 
                                       verbose=1, seed=1, nu=5/2)
  freq<-init$freq
  coef <- init$coef
  X <- expand.grid(seq(0, 1,, 101), seq(0, 1,, 101))
  colnames(X) <- c("x1", "x2")
  epsilon <- rnorm(order_tot*2)
  y <- evaluate_basis_function_FF(freq, coef, X*4)
  dfplot<- data.frame(X)
  dfplot$GP <- c(y%*% epsilon)
  library(ggplot2)
  library(viridis)
  show(ggplot(dfplot, aes(x=x1, y=x2, fill=GP))+
         geom_raster()+
         theme_bw()+
         scale_fill_viridis(option="viridis")+
         ggtitle("RFF"))
  init <- initialize_basis_function_FF(type="space-filling FF", order_tot=order_tot, dimension=2, 
                                       verbose=1, seed=1, nu=5/2)
  freq<-init$freq
  coef <- init$coef
  y <- evaluate_basis_function_FF(freq, coef, X*4)
  epsilon <- rnorm(order_tot*2)
  dfplot$GP <- c(y%*% epsilon)
  show(ggplot(dfplot, aes(x=x1, y=x2, fill=GP))+
         geom_raster()+
         theme_bw()+
         scale_fill_viridis(option="viridis")+
         ggtitle("Space-filling FF"))
  init <- initialize_basis_function_FF(type="regular", order_t=5, order_x=5, dimension=2, 
                                       verbose=1, seed=1, nu=5/2)
  freq<-init$freq
  coef <- init$coef
  y <- evaluate_basis_function_FF(freq, coef, X)
  epsilon <- rnorm(ncol(y))
  dfplot$GP <- c(y%*% epsilon)
  ggplot(dfplot, aes(x=x1, y=x2, fill=GP))+
    geom_raster()+
    theme_bw()+
    scale_fill_viridis(option="viridis")+
    ggtitle("Space-filling FF")
}

evaluate_basis_function_inducing_pt <- function(X, inducing_points,
                                                sqrtInv_covMat, kernel){
  kxX <- covMat(kernel, X=X, Xnew=inducing_points)
  functions <- kxX %*% t(sqrtInv_covMat)
  # GP <- functions %*% epsilon
  return(functions)
}

if(FALSE){
  set.seed(2)
  inducing_points <- lhsDesign(n=10, dimension=2, seed=1)$design
  colnames(inducing_points) <- c("x1", "x2")
  kernel <- kMatern(d=2, nu="5/2")
  coef(kernel) <- c(0.2, 0.2,  1)
  inputNames(kernel) <- c("x1", "x2")
  sqrtandinv_covMat <- initialize_basis_function_inducing_pt(inducing_points, kernel)
  sqrtInv_covMat <- sqrtandinv_covMat$sqrtInv_covMat
  sqrt_covMat <- sqrtandinv_covMat$sqrt_covMat
  X <- expand.grid(seq(0, 1,, 101), seq(0, 1,, 101))
  colnames(X) <- c("x1", "x2")
  epsilon <- runif(10, -2, 2)
  y <- evaluate_basis_function_inducing_pt(X, inducing_points, epsilon, sqrtInv_covMat, kernel)
  dfplot<- data.frame(X)
  dfplot$GP <- y
  dfplot2 <- data.frame(inducing_points)
  dfplot2$GP <- sqrt_covMat%*%epsilon
  library(ggplot2)
  library(viridis)
  ggplot(dfplot, aes(x=x1, y=x2, fill=GP))+
    geom_raster()+
    geom_point(data=dfplot2, pch=21, col="grey", size=3)+
    theme_bw()+
    scale_fill_viridis(option="viridis")
}

evaluate_basis_functions <- function(par_basis_functions,
                                     X,
                                     lengthscale=1){
  
  suitable_types <- c("inducing points",
                      "RFF", 
                      "space-filling FF",
                      "regular")
  type <- par_basis_functions$type
  if(!(type %in% suitable_types)){
    stop(paste0("Unsuitable value of `type`, please provide either a numerical value between 1 and ", 
                length(suitable_types), " or a string among:\n", 
                paste0(suitable_types, 
                       collapse = "; ")))
  }
  if(type=="inducing points"){
    basis_fun <- evaluate_basis_function_inducing_pt(X=t(t(X)/lengthscale), 
                                                     inducing_points=par_basis_functions$inducing_points,
                                                     sqrtInv_covMat =par_basis_functions$sqrtInv_covMat, 
                                                     kernel=par_basis_functions$kernel)
  }else{
    basis_fun <- evaluate_basis_function_FF(freq=par_basis_functions$freq,
                                            offset = par_basis_functions$offset,
                                            coef=par_basis_functions$coef,
                                            X=t(t(X)/lengthscale))
  }
  return(basis_fun)
}

initialize_basis_function_FF <- function(type=c("RFF", "space-filling FF", "regular"),
                                         nfreq=NULL, nfreq_t=NULL, nfreq_x=NULL,
                                         dimension=2, verbose=0, seed=NULL, nu=5/2){
  if(type=="RFF"){
    if(verbose==1){
      print("You selected random Fourier Features, in the current implementation, this is only available for Matérn (2k+1)/2 kernels.")
    }
    require(mvtnorm)
    freq <- mvtnorm::rmvt(n=nfreq, delta=rep(0, dimension), 
                 sigma=diag(dimension), df=2*nu)
    freq <- rbind(freq, freq)
    offset <- c(rep(0, nfreq), rep(-pi/2, nfreq))
    coef <- rep(1/sqrt(nfreq), 2*nfreq)
  }
  if(type=="space-filling FF"){
    print("You selected space-filling Fourier Features, in the current implementation, this is only available for Matérn (2k+1)/2 kernels.")
    require(DiceDesign)
    if(is.null(seed)){
      warning("Recall that DiceDesign overrides the local seed...")
    }
    X <- lhsDesign(nfreq, dimension, seed=seed)$design
    Xopt <- maximinSA_LHS(X, T0=10, c=0.99, it=2000)$design
    
    freq <- rosenblatt_transform(Xopt, dim_tot = dimension, nu=nu)
    freq <- rbind(freq, freq)
    offset <- c(rep(0, nfreq), rep(-pi/2, nfreq))
    coef <- rep(1/sqrt(nfreq), 2*nfreq)
    
  }
  if(type=="regular"){
    dim_x <- dimension - 1
    # Init
    ai_temp <- as.matrix(rbind(expand.grid(seq(1, nfreq_t),
                                           seq(0, nfreq_x)),
                               expand.grid(seq(1, nfreq_t),
                                           -seq(1, nfreq_x))))
    if(dimension > 2){
      for(i in seq(dim_x-1)){
        comb <- as.matrix(expand.grid(seq(nrow(ai_temp)),
                                      seq(-nfreq_x, nfreq_x)))
        ai_temp <- cbind(ai_temp[comb[, 1], ], comb[, 2])
      }
    }
    for(i in seq(dimension, 1)){
      ai_temp <- ai_temp[order(ai_temp[, i]), ]
    }
    ai_temp <- ai_temp*pi*2
    colnames(ai_temp)<- c("t", paste0("x", seq(dim_x)))
    nfreq <- 2*nrow(ai_temp)
    freq <- rbind(ai_temp, ai_temp)
    offset <- c(rep(0, nfreq/2), rep(-pi/2, nfreq/2))
    coef <- rep(1/sqrt(nfreq), nfreq)
  }
  return(list(freq=freq, offset=offset, coef=coef))
}

initialize_basis_function_inducing_pt <- function(inducing_points, 
                                                  # nu=c("1/2", "3/2", "5/2"),
                                                  kernel=NULL){
  require(kergp)
  dimension <- ncol(inducing_points)
  # kernel <- kMatern(d=dimension, nu=nu)
  # colnames(inducing_points) <- paste0("x", seq(dimension))
  # inputNames(kernel) <- colnames(inducing_points)
  K <- covMat(kernel, inducing_points)
  eig <- eigen(K)
  eig$values <- abs(eig$values)/2+eig$values/2
  if(any(eig$values==0)){
    stop("The inducing points/kernel provided make the covariance matrix ill-conditioned. Select different points or regularise, please.")
  }
  Ksqrt <- eig$vectors %*% diag(sqrt(eig$values)) %*% t(eig$vectors)
  Kinvsqrt <- eig$vectors %*% diag(1/sqrt(eig$values)) %*% t(eig$vectors)
  
  # for Y ~ N(0, K), Y %*% Kcholinv ~ N(0, I)
  return(list(sqrtInv_covMat=Kinvsqrt, sqrt_covMat=Ksqrt))
}

initialize_basis_functions <- function(par_basis_functions,
                                       verbose=0,
                                       seed=NULL){
  type <- par_basis_functions$type
  
  suitable_types <- c("inducing points",
                      "RFF", 
                      "space-filling FF",
                      "regular")
  if(type %in% seq(4)){
    warning(paste0("You provided `type` as a numeric variable with value: ", 
                   type, 
                   ". It will be interpreted as: `", suitable_types[type], "`."))
    type <- suitable_types[type]
    par_basis_functions$type <<- type
  }
  if(!(type %in% suitable_types)){
    stop(paste0("Unsuitable value of `type`, please provide either a numerical value between 1 and ", 
                length(suitable_types), " or a string among:\n", 
                paste0(suitable_types, 
                       collapse = "; ")))
  }
  if(type=="inducing points"){
    init <- initialize_basis_function_inducing_pt(par_basis_functions$inducing_points, 
                                                  kernel=par_basis_functions$kernel)
    par_basis_functions$sqrt_covMat <<- init$sqrt_covMat
    par_basis_functions$sqrtInv_covMat <<- init$sqrtInv_covMat
    par_basis_functions$order_tot <<- nrow(inducing_points)
    par_basis_functions$dimension <<- ncol(inducing_points)
  }else{
    if(type=="regular"){
      init <- initialize_basis_function_FF(type=type,
                                           nfreq_t=par_basis_functions$nfreq_t,
                                           nfreq_x=par_basis_functions$nfreq_x,
                                           dimension = par_basis_functions$dimension,
                                           verbose=verbose,
                                           seed=seed)
      par_basis_functions$freq <<- init$freq
      par_basis_functions$coef <<- init$coef
      par_basis_functions$offset <<- init$offset
      par_basis_functions$order_tot <<- length(init$coef)
    }else{
      init <- initialize_basis_function_FF(type=type,
                                           nfreq=par_basis_functions$nfreq,
                                           dimension = par_basis_functions$dimension,
                                           verbose=verbose,
                                           seed=seed,
                                           nu=par_basis_functions$nu)
      
      par_basis_functions$freq <<- init$freq
      par_basis_functions$coef <<- init$coef
      par_basis_functions$offset <<- init$offset
      par_basis_functions$order_tot <<- length(init$coef)
    }
  }
}

if(FALSE){
  library(kergp)
  set.seed(2)
  inducing_points <- lhsDesign(n=10, dimension=2, seed=1)$design
  colnames(inducing_points) <- c("x1", "x2")
  kernel <- kMatern(d=2, nu="5/2")
  coef(kernel) <- c(0.2, 0.2,  1)
  inputNames(kernel) <- c("x1", "x2")
  
  par_basis_functions <- list(type="inducing points",
                              inducing_points=inducing_points,
                              kernel=kernel)
  summary(par_basis_functions)
  initialize_basis_functions(par_basis_functions, verbose=1)
  
  summary(par_basis_functions)
  
  
  
  par_basis_functions <- list(type="regular",
                              dimension=2,
                              nfreq_t =3, 
                              nfreq_x =3)
  
  initialize_basis_functions(par_basis_functions, verbose=1)
  par_basis_functions <- list(type="RFF",
                              dimension=2,
                              nfreq = 20,
                              nu=5/2)
  initialize_basis_functions(par_basis_functions, verbose=1)
}

heuristic_find_variance <- function(par_basis_functions, nsimu=1000, lengthscale=1,
                                    grid_size=101, plot=FALSE){
  dimension <- par_basis_functions$dimension
  order_tot <- par_basis_functions$order_tot
  type <- par_basis_functions$type
  epsilon <- matrix(rnorm(nsimu*order_tot), ncol=order_tot)
  X <- expand.grid(rep(list(seq(0, 1,, grid_size)), dimension))
  X <- t(t(X)/lengthscale)
  if(par_basis_functions$type=="inducing points"){
    colnames(X)<- inputNames(par_basis_functions$kernel)
  }
  funX <- evaluate_basis_functions(par_basis_functions, X)
  GPX <- funX%*% t(epsilon)
  range_size <- apply(GPX, 2, function(x){
    return(diff(range(x)))
  })
  if(plot){
    require(ggplot2)
    show(ggplot(data.frame(x=range_size), aes(x=x))+
           geom_histogram(mapping=aes(y=..density..),bins=30, alpha=0.1, col="black")+
           geom_density(lty=2)+
           geom_vline(xintercept = mean(range_size), col="blue", lwd=2)+
           geom_vline(xintercept = median(range_size), col="darkgreen", lwd=2)+
           theme_bw()+
           labs(caption=paste0("Mean: ", round(mean(range_size), 2), " (blue).\n",
                               "Median: ", round(mean(range_size), 2), " (green)."))+
           ggtitle("Range of values (max-min) of the GP"))
  }
  return(range_size)
}

