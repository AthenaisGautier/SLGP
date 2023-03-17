rosenblatt_transform <- function(x, dim_tot, nu=5/2){
  ## Starts from a x presumably uniform in [0, 1]. x is either a n*dim_tot matrix or a dim_tot-vector
  if(is.vector(x)){
    x <- matrix(x, nrow=1)
  }
  if(!is.matrix(x)){
    stop("Please, apply the Rosenblatt transform to a `x` that is either a vector or a matrix")
  }
  if(ncol(x)!=dim_tot){
    stop("Check `x`'s dimension for the Rosenblatt tranform.")
  }
  # The first marginal is computed:
  newx <- 0*x
  newx[, 1] <- qt(x[, 1], df=2*nu)
  partial_sum <- 0
  for(i in seq(2, dim_tot)){
    #Conditionals
    partial_sum <- partial_sum + x[, i-1]^2
    newx[, i] <- qt(x[, i], df=2*nu+i-1)
    newx[, i] <- newx[, i]*sqrt(2*nu+partial_sum)/sqrt(2*nu +i -1)
  }
  return(newx)
}

if(FALSE){
  n <- 20
  seed <- 3
  set.seed(seed)
  library(ggplot2)
  library(mvtnorm)
  library(viridis)
  library(tidyverse)
  u <- seq(-4, 4, 0.05)
  dim_tot <- 2
  nu <- 5/2
  data <- data.frame(expand.grid(u, u))
  colnames(data) <- c("w1", "w2")
  data$density <- mvtnorm::dmvt(x=data, delta=rep(0, dim_tot),  
                       sigma=diag(dim_tot), df=2*nu, log=FALSE)
  colmap <- c(NA, rev(viridis(n=10, option="magma")))
  plot1<- ggplot(mapping=aes(x=w1, y=w2))+
    geom_raster(data, mapping=aes(fill=density))+
    theme_bw()+
    scale_fill_gradientn(
      name = "Probability density",
      colors=colmap,
      # here we use guide_colourbar because it is still a continuous scale
      guide = guide_colorbar(
        direction = "horizontal",
        barheight = unit(2, units = "mm"),
        barwidth = unit(50, units = "mm"),
        draw.ulim = F,
        title.position = 'top',
        # some shifting around
        title.hjust = 0.5,
        label.hjust = 0.5
      ))+
    theme(legend.position = "none")
  samples <- data.frame(mvtnorm::rmvt(n=n, delta=rep(0, dim_tot),  
                  sigma=diag(dim_tot), df=2*nu))
  colnames(samples) <- c("w1", "w2")
  samples$method <- "Random"
  library(DiceDesign)
  X <- lhsDesign(n, dim_tot, seed=seed)$design
  Xopt <- maximinSA_LHS(X, T0=10, c=0.99, it=2000)$design
  
  newsamp <- data.frame(rosenblatt_transform(Xopt,
                                  dim_tot = dim_tot, nu=nu))
  colnames(newsamp) <- c("w1", "w2")
  newsamp$method <- "Space-filling"
  samples <- rbind(samples, newsamp)
  plot1 +
    geom_point(data=samples, pch=21, col="black", fill="grey")+
    facet_grid(.~method)+
    xlab(expression(omega[1]))+
    ylab(expression(omega[2]))
  ggsave("./Figures/illustration_space_fill_RFF_20.png", width=6, height=3)
  
  set.seed(seed+4)
  df_plot <- data.frame(expand.grid(seq(0, 5,, 101), seq(0, 5,, 101)))
  colnames(df_plot) <- c("x", "t")
  temp <- as.matrix(df_plot) %*% t(samples[, 1:2])
  epsilon <- rnorm(2*n)
  temp <- cbind(cos(temp), sin(temp))
  temp <- temp[, c(1:n, (1:n)+2*n, (1:n)+n, (1:n)+3*n)]
  df_plot <- rbind(df_plot, df_plot)
  df_plot$simu1 <-  c(temp[, 1:(2*n)] %*% epsilon, 
                      temp[, 2*n+1:(2*n)] %*% epsilon)/sqrt(2*n)
  epsilon <- rnorm(2*n)
  df_plot$simu2 <-  c(temp[, 1:(2*n)] %*% epsilon, 
                      temp[, 2*n+1:(2*n)] %*% epsilon)/sqrt(2*n)
  epsilon <- rnorm(2*n)
  df_plot$simu3 <-  c(temp[, 1:(2*n)] %*% epsilon, 
                      temp[, 2*n+(1:(2*n))] %*% epsilon)/sqrt(2*n)
  colnames(df_plot) <- c("x", "t", "Simulated SLGP: 1",
                         "Simulated SLGP: 2", "Simulated SLGP: 3")
  df_plot$method <- c(rep("Random", nrow(df_plot)/2), 
                      rep("Space-Filling", nrow(df_plot)/2))
  df_plot%>%
    pivot_longer(-c("x", "t", "method"))%>%
    group_by(x, method, name)%>%
    mutate(value=exp(value))%>%
    mutate(value=value/mean(value)/5)%>%
    ggplot(aes(x=x, y=t, fill=value))+
    geom_raster()+
    facet_grid(method~name)+
    scale_fill_gradientn(
      name = "Probability density",
      colors=colmap,
      # here we use guide_colourbar because it is still a continuous scale
      guide = guide_colorbar(
        direction = "horizontal",
        barheight = unit(2, units = "mm"),
        barwidth = unit(50, units = "mm"),
        draw.ulim = F,
        title.position = 'top',
        # some shifting around
        title.hjust = 0.5,
        label.hjust = 0.5
      ))+
    theme_bw()+
    theme(legend.position = "bottom")
}
require(stringr)
give_names_columns <- function(range_min, range_max, charac){
  return(paste0(charac,
                str_sub(paste0("00000", range_min:range_max), start=-5)))
}
print_time <- function(time){
  seconds <- as.numeric(time, units="secs")
  min <- seconds %/% 60
  seconds <- seconds - 60*min
  hours <- min %/% 60
  min <- min - 60*hours
  if(hours==0){
    chr_h <- ""
  }else{
    chr_h <- paste0(hours, "h")
  }
  if(min==0){
    chr_m <- ""
  }else{
    chr_m <- paste0(min, "m")
  }
  chr_s <- paste0(round(seconds), "s")
  return(paste0(chr_h, chr_m, chr_s))
}

