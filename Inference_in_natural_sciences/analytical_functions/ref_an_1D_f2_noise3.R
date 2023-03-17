library(triangle)

ref_field <- function(X){
  res <- data.frame(X)
  res <- res[, c("t", "x1")]
  med <- ref_median(x=res[, c(2)], scaled=T)
  alpha <- 0.8-50*res[, c(2)]*(res[, c(2)]-0.7)^4
  b1 <- med
  b2 <- (1-med)
  b3 <- pmin(b2, b1)
  res$pdf <- sapply(seq(length(med)), FUN=function(i){
    alpha[i]*(0.5*dtriangle(X$t[i]-med[i], -1*b3[i], 0*b3[i], -0.5*b3[i])+
                0.5*dtriangle(X$t[i]-med[i], 0*b3[i], 1*b3[i], 0.5*b3[i])) +
      (1-alpha[i])*dtriangle(X$t[i]-med[i], -0.5*b3[i], 0.5*b3[i], 0*b3[i])})
  res$cdf <- sapply(seq(length(med)), FUN=function(i){
    alpha[i]*(0.5*ptriangle(X$t[i]-med[i], -1*b3[i], 0*b3[i], -0.5*b3[i])+
                0.5*ptriangle(X$t[i]-med[i], 0*b3[i], 1*b3[i], 0.5*b3[i])) +
      (1-alpha[i])*ptriangle(X$t[i]-med[i], -0.5*b3[i], 0.5*b3[i], 0*b3[i])})
  return(res)
}

ref_median <- function(x, scaled=T){
  x <- c(x)
  y <- x*10-5
  med <- (y^2-5*y+6)/(y^2+1)
  med <- (med+0.1)/(7.1+0.1)
  med <- med*0.7 + 0.15
  return(med)
}

sampling_values <- function(unique_loc, unique_freq){
  unique_loc <- c(unique_loc)
  med <- ref_median(x=unique_loc, scaled=T)
  alpha <- 0.8-50*unique_loc*(unique_loc-0.7)^4
  
  b1 <- med
  b2 <- (1-med)
  b3 <- pmin(b2, b1)
  
  unique_freq <- c(unique_freq)
  n <- length(unique_freq)
  samples <- data.frame(t=1, x1=1)[0, ]
  for(i in seq(n)){
    r1 <- runif(unique_freq[i])
    r2 <- runif(unique_freq[i])
    new_samp <- data.frame(t = med[i]+
                             (r1 <= alpha[i])*((r2 <= 0.5) * rtriangle(unique_freq[i], -1*b3[i], 0*b3[i], -0.5*b3[i]) +
                                                 (r2 > 0.5) * rtriangle(unique_freq[i], 0*b3[i], 1*b3[i], 0.5*b3[i])) +                           +
                             (r1 > alpha[i])*rtriangle(unique_freq[i], -0.5*b3[i], 0.5*b3[i], 0*b3[i]),
                           x1 = unique_loc[i])
    samples <- rbind(samples, new_samp)
  }
  return(samples)
}

