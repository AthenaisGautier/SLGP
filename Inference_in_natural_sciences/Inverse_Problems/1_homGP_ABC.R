# setwd("C:/Users/athen/Research/Codes/repo_PHD/Inference_in_natural_sciences/Inverse_Problems/")
setwd("/storage/homefs/ag19v929/Research/Codes/repo_PHD/Inference_in_natural_sciences/Inverse_Problems")

args <- commandArgs(trailingOnly=TRUE)
job_number <-  as.numeric(args[1])
# job_number <- 1
sampling <- as.character(args[2])
# sampling <- "RandomUnif"
# sampling <- "RandomBatch" 
# sampling <- "EIV" 
# sampling <- "EIMAD" 
config <- read.csv("config.txt")
config <- config[job_number, ]

library(sn)
library(conflicted)
conflict_prefer("sd", "stats") 
conflict_prefer("mle", "kergp") 
conflict_prefer("dst", "LaplacesDemon") 
conflict_prefer("pst", "LaplacesDemon") 
conflict_prefer("qst", "LaplacesDemon") 
conflict_prefer("rst", "LaplacesDemon") 
conflict_prefer("tr", "LaplacesDemon") 

set.seed(config$rep)

if(config$case=="geol"){
  load(file=paste0("../ref_field_geol_", config$geol, "_source_", config$source, ".RData"))
}else{
  source(paste0("../analytical_functions/ref_an_1D_f", config$fmed, "_noise", config$noise, ".R"))
  df_ref <- data.frame(expand.grid(seq(0, 1,, 101), seq(0, 1,, 101)))
  colnames(df_ref) <- c("t", "x1")
  df_ref <- ref_field(df_ref)
}
treshold <- 0.15
df_ref <- round(df_ref, 14)
ref_ABC <- df_ref[df_ref$t==treshold, c("x1", "cdf")]
ref_ABC$cdf <- ref_ABC$cdf/mean(ref_ABC$cdf)
colnames(ref_ABC)[2] <- "ABCposterior"
if(config$case=="geol"){
  samples_geol_ABC <- samples_geol[samples_geol$t<= treshold, ]
}

candidates <- seq(0, 1,, 51)
nmax <- 1500
nstart <- 50
loc_start <- sample(x=candidates, nstart, replace=TRUE)
r_start <- rank(loc_start, ties.method = "first")
loc_start <- sort(loc_start)
source("../utils.R")
tab_loc <- as.data.frame(table(loc_start))
tab_loc$loc_start <- as.numeric(as.character(tab_loc$loc_start))
if(config$case=="geol"){
  sample_start <- data.frame(t= unlist(sapply(seq(nrow(tab_loc)), function(i){
    sub_df <- df_ref[round(tab_loc[i, 1], 14)== df_ref$x1, ]
    return(sample_from_custom_cdf(n=tab_loc[i, 2], t=sub_df$t, cdf=sub_df$cdf))
  })),
  x1=loc_start)
}else{
  source(paste0("../analytical_functions/ref_an_1D_f", config$fmed, "_noise", config$noise, ".R"))
  sample_start <-  data.frame(t= unlist(sapply(seq(nrow(tab_loc)), function(i){
    temp <- sampling_values(unique_loc=tab_loc[i, 1], unique_freq=tab_loc[i, 2])
    return(temp$t)
  })),
  x1=loc_start)
}
sample_start <- sample_start[r_start, ]

samples <- sample_start

library(kergp)
cov <- kMatern(d=1, nu="5/2")
coef(cov) <- c(0.05, 1)
inputNames(cov) <- c("x1")
mod <- gp(t~., data=samples, cov=cov, 
          noise = TRUE,
          varNoiseLower = 1e-5,
          varNoiseIni=0.01,
          varNoiseUpper=0.1)
noiIni <- mod$varNoise

df_ABC <- data.frame(x1=seq(0, 1,, 101))

pred <- predict(mod, df_ABC, covCompute = TRUE, seCompute = TRUE, type="UK")
library(MASS)
n_real <- 100
simu_post_GP <- t(mvrnorm(n=n_real, mu=pred$mean, Sigma=pred$cov))
colnames(simu_post_GP) <- paste0("simu_", seq(n_real))
df_ABC <- cbind(df_ABC, pnorm(treshold-simu_post_GP, sd=sqrt(mod$varNoise)))
for(i in seq(n_real)){
  df_ABC[, i+1] <- df_ABC[, i+1] / mean(df_ABC[, i+1])
}
df_ABC$mean_ABCposterior <- as.vector(pnorm(treshold, mean=pred$mean, sd=pred$sd))
df_ABC$mean_ABCposterior <- df_ABC$mean_ABCposterior/mean(df_ABC$mean_ABCposterior)

K_QQ <- K_QP <- K_PP <- exp(-as.matrix(dist(df_ABC$x1, upper=TRUE, diag=TRUE))^2/0.1)
values_Q <- df_ABC$mean_ABCposterior
values_P <- ref_ABC$ABCposterior
score <- score_pdf(values_Q=values_Q, values_P=values_P,
                   K_PP=K_PP, K_QP=K_QP, K_QQ=K_QQ, 
                   discrete_Q =FALSE, discrete_P = FALSE)
# The probabilistic pred
values_Q <- df_ABC[, colnames(simu_post_GP)]
score2 <- score_pdf(values_Q=values_Q, values_P=values_P, 
                    K_PP=K_PP, K_QP=K_QP, K_QQ=K_QQ, 
                    discrete_Q =FALSE, discrete_P = FALSE)

# plot(df_ABC, type="l")
score_val <- c(score, score2)

if(config$case=="geol"){
  K_PP2 <- exp(-as.matrix(dist(samples_geol_ABC$x1, upper=TRUE, diag=TRUE))^2/0.1)
  K_QP2 <- outer(df_ABC$x1, samples_geol_ABC$x1,
                 Vectorize(function(x, y){return(exp(-sqrt(sum((x-y)^2))/0.1))}))
  values_P2 <- rep(1, nrow(samples_geol_ABC))
  
  # Just the mean pred
  values_Q <- df_ABC$mean_ABCposterior
  score <- score_pdf(values_Q=values_Q, values_P=values_P2,
                     K_PP=K_PP2, K_QP=K_QP2, K_QQ=K_QQ, 
                     discrete_Q =FALSE, discrete_P = TRUE)
  # The probabilistic pred
  values_Q <- df_ABC[, colnames(simu_post_GP)]
  score2 <- score_pdf(values_Q=values_Q, values_P=values_P2, 
                      K_PP=K_PP2, K_QP=K_QP2, K_QQ=K_QQ, 
                      discrete_Q =FALSE, discrete_P = TRUE)
  score_val2 <- c(score, score2)
}

if(FALSE){
  library(ggplot2)
  library(gridExtra)
  plot1 <- ggplot()+
    geom_point(data=samples, aes(x=x1, y=t))+
    geom_line(data=data.frame(x1=seq(0, 1,, 101), y=pred$mean), aes(x=x1, y=y))+
    geom_ribbon(data=data.frame(x1=seq(0, 1,, 101), ymin=pred$lower95, ymax=pred$upper95),
                aes(x=x1, ymin=ymin, ymax=ymax), alpha=0.2, col="blue", lty=2)+
    theme_bw()
  plot2 <- ggplot(data.frame(n=nstart+sizebatch*(0:nbatch)[1:length(score_val)],
                             score=score_val),
                  aes(x=n, y=score))+
    geom_line()+
    theme_bw()
  grid.arrange(plot1, plot2, nrow=2, heights=c(4, 2))
}

sizebatch <- 25
nbatch <- (nmax-nstart)/sizebatch
res <- data.frame(config, "PDF VS mean", t(c(score_val[1], rep(NA, nbatch))[1:(nbatch+1)]))
temp <- data.frame(config, 
                   "PDF VS prob", 
                   t(c(score_val[2], rep(NA, nbatch))[1:(nbatch+1)]))

colnames(res)[-c(1:ncol(config))]<- c("reference", paste0("Score_", nstart+sizebatch*(0:nbatch)))
colnames(temp)<- colnames(res)
res <-rbind(res,temp)
if(config$case=="geol"){
  temp <-data.frame(config, "sample VS mean", t(c(score_val2[1], rep(NA, nbatch))[1:(nbatch+1)]))
  colnames(temp)<- colnames(res)
  res <-rbind(res,temp)
  temp <-data.frame(config, "sample VS prob", t(c(score_val2[2], rep(NA, nbatch))[1:(nbatch+1)]))
  colnames(temp)<- colnames(res)
  res <-rbind(res,temp)
}

ibatch<-1
for(ibatch in seq(nbatch)){
  set.seed(config$rep*+1000*ibatch)
  print(nstart+ibatch*sizebatch)
  if(sampling=="RandomUnif"){
    loc_new <- sample(x=candidates, sizebatch, replace=TRUE)
    r_new <- rank(loc_new, ties.method = "first")
    loc_new <- sort(loc_new)
    source("../utils.R")
    tab_loc <- as.data.frame(table(loc_new))
    tab_loc$loc_new <- as.numeric(as.character(tab_loc$loc_new))
    if(config$case=="geol"){
      sample_new <- data.frame(t= unlist(sapply(seq(nrow(tab_loc)), function(i){
        sub_df <- df_ref[round(tab_loc[i, 1], 14)== df_ref$x1, ]
        return(sample_from_custom_cdf(n=tab_loc[i, 2], t=sub_df$t, cdf=sub_df$cdf))
      })),
      x1=loc_new)
    }else{
      source(paste0("../analytical_functions/ref_an_1D_f", config$fmed, "_noise", config$noise, ".R"))
      sample_new <-  data.frame(t= unlist(sapply(seq(nrow(tab_loc)), function(i){
        temp <- sampling_values(unique_loc=tab_loc[i, 1], unique_freq=tab_loc[i, 2])
        return(temp$t)
      })),
      x1=loc_new)
    }
  }
  if(sampling=="RandomBatch"){
    loc_new <- sample(x=candidates, 1)
    source("../utils.R")
    if(config$case=="geol"){
      sub_df <- df_ref[round(loc_new, 14)== df_ref$x1, ]
      sample_new <- data.frame(t= sample_from_custom_cdf(n=sizebatch, t=sub_df$t, cdf=sub_df$cdf),
                               x1=loc_new)
      rm(sub_df)
    }else{
      source(paste0("../analytical_functions/ref_an_1D_f", config$fmed, "_noise", config$noise, ".R"))
      sample_new <-  sampling_values(unique_loc=loc_new, unique_freq=sizebatch)
    }
  }
  if(sampling=="EIV"){
    an <- as.vector((treshold-pred$mean)/(pred$sd^2+mod$varNoise))
    crit <- -sapply(candidates, function(x){
      indx <- which.min(abs(df_ABC$x1-x))
      vect <- pred$cov[rep(indx, sizebatch), ]
      cn <- diag(t(vect)%*% solve(pred$sd[indx]^2+diag(mod$varNoise, sizebatch)) %*% vect)
      T1 <- sapply(seq(length(an)), function(i){
        T.Owen(an[i], sqrt(mod$varNoise+ pred$sd[i]^2-cn[i])/ 
                 sqrt(mod$varNoise+ pred$sd[i]^2+cn[i]))
      })
      T2 <- sapply(seq(length(an)), function(i){
        T.Owen(an[i], sqrt(mod$varNoise)/ sqrt(mod$varNoise+ 2*pred$sd[i]^2))
      })
      return(2*mean(T1-T2))
    })
    loc_new <- candidates[which.min(crit)]
    source("../utils.R")
    if(config$case=="geol"){
      sub_df <- df_ref[round(loc_new, 14)== df_ref$x1, ]
      sample_new <- data.frame(t= sample_from_custom_cdf(n=sizebatch, t=sub_df$t, cdf=sub_df$cdf),
                               x1=loc_new)
      rm(sub_df)
    }else{
      source(paste0("../analytical_functions/ref_an_1D_f", config$fmed, "_noise", config$noise, ".R"))
      sample_new <-  sampling_values(unique_loc=loc_new, unique_freq=sizebatch)
    }
  }
  if(sampling=="EIMAD"){
    an <- as.vector((treshold-pred$mean)/(pred$sd^2+mod$varNoise))
    crit <- -sapply(candidates, function(x){
      indx <- which.min(abs(df_ABC$x1-x))
      vect <- pred$cov[rep(indx, sizebatch), ]
      cn <- diag(t(vect)%*% solve(pred$sd[indx]^2+diag(mod$varNoise, sizebatch)) %*% vect)
      T1 <- sapply(seq(length(an)), function(i){
        T.Owen(an[i], sqrt(pred$sd[i]^2-cn[i])/ 
                 sqrt(pred$sd[i]^2+cn[i]))
      })
      
      return(2*mean(T1))
    })
    loc_new <- candidates[which.min(crit)]
    
    source("../utils.R")
    if(config$case=="geol"){
      sub_df <- df_ref[round(loc_new, 14)== df_ref$x1, ]
      sample_new <- data.frame(t= sample_from_custom_cdf(n=sizebatch, t=sub_df$t, cdf=sub_df$cdf),
                               x1=loc_new)
      rm(sub_df)
    }else{
      source(paste0("../analytical_functions/ref_an_1D_f", config$fmed, "_noise", config$noise, ".R"))
      sample_new <-  sampling_values(unique_loc=loc_new, unique_freq=sizebatch)
    }
  }
  # show(plot(candidates, crit, "l"))
  samples <- rbind(samples, sample_new)
  
  mod <- gp(t~., data=samples, cov=cov, 
            noise = TRUE,
            varNoiseLower = 1e-5,
            varNoiseIni=noiIni,
            varNoiseUpper=0.1)
  noiIni <- mod$varNoise
  
  
  df_ABC <- data.frame(x1=seq(0, 1,, 101))
  pred <- predict(mod, df_ABC, covCompute = TRUE, seCompute = TRUE, type="UK")
  colnames(simu_post_GP) <- paste0("simu_", seq(n_real))
  df_ABC <- cbind(df_ABC, pnorm(treshold-simu_post_GP, sd=sqrt(mod$varNoise)))
  for(i in seq(n_real)){
    df_ABC[, i+1] <- df_ABC[, i+1] / mean(df_ABC[, i+1])
  }
  df_ABC$mean_ABCposterior <- as.vector(pnorm(treshold, mean=pred$mean, sd=pred$sd))
  df_ABC$mean_ABCposterior <- df_ABC$mean_ABCposterior/mean(df_ABC$mean_ABCposterior)
  values_Q <- df_ABC$mean_ABCposterior
  score <- score_pdf(values_Q=values_Q, values_P=values_P,
                     K_PP=K_PP, K_QP=K_QP, K_QQ=K_QQ, 
                     discrete_Q =FALSE, discrete_P = FALSE)
  # The probabilistic pred
  values_Q <- df_ABC[, colnames(simu_post_GP)]
  score2 <- score_pdf(values_Q=values_Q, values_P=values_P, 
                      K_PP=K_PP, K_QP=K_QP, K_QQ=K_QQ, 
                      discrete_Q =FALSE, discrete_P = FALSE)
  
  # plot(df_ABC, type="l")
  score_val <- c(score_val, score, score2)
  res[, ncol(config)+ibatch+2] <- c(score, score2)
  
  if(config$case=="geol"){
    # Just the mean pred
    values_Q <- df_ABC$mean_ABCposterior
    score <- score_pdf(values_Q=values_Q, values_P=values_P2,
                       K_PP=K_PP2, K_QP=K_QP2, K_QQ=K_QQ, 
                       discrete_Q =FALSE, discrete_P = TRUE)
    # The probabilistic pred
    values_Q <- df_ABC[, colnames(simu_post_GP)]
    score2 <- score_pdf(values_Q=values_Q, values_P=values_P2, 
                        K_PP=K_PP2, K_QP=K_QP2, K_QQ=K_QQ, 
                        discrete_Q =FALSE, discrete_P = TRUE)
    score_val2 <- c(score_val2, score, score2)
    res[3:4, ncol(config)+ibatch+2] <- c(score, score2)
  }
  
  if(FALSE){
    library(ggplot2)
    library(gridExtra)
    plot1 <- ggplot()+
      geom_point(data=samples, aes(x=x1, y=t))+
      geom_line(data=data.frame(x1=seq(0, 1,, 101), y=pred$mean), aes(x=x1, y=y))+
      geom_ribbon(data=data.frame(x1=seq(0, 1,, 101), ymin=pred$lower95, ymax=pred$upper95),
                  aes(x=x1, ymin=ymin, ymax=ymax), alpha=0.2, col="blue", lty=2)+
      theme_bw()
    plot2 <- ggplot(data.frame(n=nstart+sizebatch*(0:nbatch)[1:length(score_val)],
                               score=score_val),
                    aes(x=n, y=score))+
      geom_line()+
      theme_bw()
    grid.arrange(plot1, plot2, nrow=2, heights=c(4, 2))
  }
  
  write.csv(res,
            file=paste0("./res/scores/homGP_sampling_", sampling, "_jobnumber_", job_number, ".txt"))
  
  write.csv(samples,
            file=paste0("./res/samples/homGP_sampling_", sampling, "_jobnumber_", job_number, ".txt"))
}

