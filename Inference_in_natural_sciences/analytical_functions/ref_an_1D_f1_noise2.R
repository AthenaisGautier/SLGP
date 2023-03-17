library(truncnorm)

ref_field <- function(X){
  res <- data.frame(X)
  res <- res[, c("t", "x1")]
  med <- ref_median(x=res[, c(2)], scaled=T)
  b1 <- med
  b2 <- (1-med)
  b3 <- pmin(b2, b1)
  res$pdf <- sapply(seq(length(med)), FUN=function(i){
    dtruncnorm(X$t[i]-med[i], a=-b3[i], b=b3[i], mean=0, sd=b3[i]/3)
  })
  res$cdf <- sapply(seq(length(med)), FUN=function(i){
    ptruncnorm(X$t[i]-med[i], a=-b3[i], b=b3[i], mean=0, sd=b3[i]/3)
  })
  return(res)
}

ref_median <- function(x, scaled=T){
  x <- c(x)
  y <- x*(7.5-2.7)+2.7
  med <- sin(y)+sin(10*y/3)
  med <- (med+1.9)/(1.9+0.9)
  med <- med*0.7 + 0.15
  return(med)
}

sampling_values <- function(unique_loc, unique_freq){
  unique_loc <- c(unique_loc)
  med <- ref_median(x=unique_loc, scaled=T)
  b1 <- med
  b2 <- (1-med)
  b3 <- pmin(b2, b1)
  unique_freq <- c(unique_freq)
  n <- length(unique_freq)
  samples <- data.frame(t=1, x1=1)[0, ]
  for(i in seq(n)){
    new_samp <- data.frame(t = med[i]+rtruncnorm(unique_freq[i], 
                                                 a=-b3[i], b=b3[i], 
                                                 mean=0, sd=b3[i]/3),
                           x1 = unique_loc[i])
    samples <- rbind(samples, new_samp)
  }
  return(samples)
}

