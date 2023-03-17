util_fun_median <- function(t, cdf, q){
  f_temp <- approxfun(x=t, y=cdf)
  f_inv <- inverse(f_temp, lower = 0, upper = 1)
  return(f_inv(q))
}

library(GoFKernel)

sample_from_custom_cdf <- function(n, t, cdf){
  r <- runif(n)
  f_temp <- approxfun(x=t, y=cdf)
  f_inv <- inverse(f_temp, lower = 0, upper = 1)
  samp <- sapply(r, f_inv)
  return(samp)
}

score_pdf <- function(values_Q, values_P, K_PP, K_QP, K_QQ, discrete_Q =FALSE,
                      discrete_P = FALSE){
  if(is.null(dim(values_Q))){
    values_Q <- matrix(values_Q)
  }else{
    values_Q <- as.matrix(values_Q)
  }
  if(is.null(dim(values_P))){
    values_P <- matrix(values_P)
  }else{
    values_P <- as.matrix(values_P)
  }
  np <- ncol(values_P)
  nq <- ncol(values_Q)
  
  if(discrete_P){
    term_P <- 2*mean(lower.tri(K_PP))
  }else{
    if(np>1){
      term_P <- mean(unlist(sapply(seq(2, np), function(i){
        sapply(seq(i-1), function(j){
          mean(values_P[, i]*t(values_P[, j]*K_PP))
        })
      })))
    }else{
      term_P <- mean(values_P[,1]*t(values_P[,1]*K_PP))
    }
  }
  if(discrete_Q){
    term_Q <- 2*mean(lower.tri(K_QQ))
  }else{
    if(nq>1){
      term_Q <- mean(unlist(sapply(seq(2, nq), function(i){
        sapply(seq(i-1), function(j){
          mean(values_Q[, i]*t(values_Q[, j]*K_QQ))
        })
      })))
    }else{
      term_Q <- mean(values_Q[,1]*t(values_Q[,1]*K_QQ))
    }
  }
  term_QP <- mean(sapply(seq(nq), function(i){
    sapply(seq(np), function(j){
      mean(values_P[, j]*t(values_Q[, i]*K_QP))
    })
  }))
  return(term_P+term_Q-2*term_QP)
}
