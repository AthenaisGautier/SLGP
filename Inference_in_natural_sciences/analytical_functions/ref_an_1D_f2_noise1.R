library(truncnorm)

ref_field <- function(X){
  res <- data.frame(X)
  res <- res[, c("t", "x1")]
  med <- ref_median(x=res[, c(2)], scaled=T)
  res$pdf <- sapply(seq(length(med)), FUN=function(i){
    dtruncnorm(X$t[i]-med[i], a=-0.15, b=0.15, mean=0, sd=0.05)
  })
  res$cdf <- sapply(seq(length(med)), FUN=function(i){
    ptruncnorm(X$t[i]-med[i], a=-0.15, b=0.15, mean=0, sd=0.05)
  })
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
  unique_freq <- c(unique_freq)
  n <- length(unique_freq)
  samples <- data.frame(t=1, x1=1)[0, ]
  for(i in seq(n)){
    new_samp <- data.frame(t = med[i]+rtruncnorm(unique_freq[i], a=-0.15, b=0.15, mean=0, sd=0.05),
                           x1 = unique_loc[i])
    samples <- rbind(samples, new_samp)
  }
  return(samples)
}

# # Real min
# u <- seq(0, 1,, 101)
# y <- ref_median(u)
# range(y)
# opt <- optim(fn=ref_median, lower=0, upper=1, method="L-BFGS-B", par=u[which.min(y)])
# D <- 1
# fun <- 2
# write.csv(file=paste0("opt_an_", D, "D_f", fun, ".txt"),
#           data.frame(x=c(opt$par[1]), value=opt$value))
