# setwd("C:/Users/athen/Research/Codes/repo_PHD/Inference_in_natural_sciences/Inverse_Problems/")
setwd("/storage/homefs/ag19v929/Research/Codes/repo_PHD/Inference_in_natural_sciences/Inverse_Problems")
args <- commandArgs(trailingOnly=TRUE)
job_number <-  as.numeric(args[1])
# job_number <- 1
sampling <- as.character(args[2])
# sampling <- "RandomUnif"
config <- read.csv("config.txt")
config <- config[job_number, ]


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
nmax <- 10000
nstart <- 50
loc_start <- runif(nstart)
r_start <- rank(loc_start, ties.method = "first")
loc_start <- sort(loc_start)
source("../utils.R")
tab_loc <- as.data.frame(table(loc_start))
tab_loc$loc_start <- as.numeric(as.character(tab_loc$loc_start))
library(pracma)
u <- unique(df_ref$t)
if(config$case=="geol"){
  sample_start <- data.frame(t= unlist(sapply(seq(nrow(tab_loc)), function(i){
    sub_df <- df_ref[abs(round(tab_loc[i, 1], 14)- df_ref$x1)<=1/100, ]
    value_cdf <- interp2(x=unique(sub_df$x1),
            y=u,
            Z=matrix(sub_df$cdf, ncol=2),
            yp=u, 
            xp=rep(tab_loc[i, 1], length(u)))
    return(sample_from_custom_cdf(n=tab_loc[i, 2], t=u, cdf=value_cdf))
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
samples <- sample_start[r_start, ]

sample_ABC <- samples[samples$t <= treshold, ]

values_P <- ref_ABC$ABCposterior
K_PP <- exp(-as.matrix(dist(ref_ABC$x1, upper=TRUE, diag=TRUE))^2/0.1)

if(config$case=="geol"){
  #Compare against sample as well
  K_PP2 <- exp(-as.matrix(dist(samples_geol_ABC$x1, upper=TRUE, diag=TRUE))^2/0.1)
  values_P2 <- rep(1, nrow(samples_geol_ABC))
}
if(nrow(sample_ABC)==0){
  score <- NA
  if(config$case=="geol"){
    score2 <- NA
    score_val2 <- c(score2)
  }
}else{
  K_QQ <- exp(-as.matrix(dist(sample_ABC$x1, upper=TRUE, diag=TRUE))^2/0.1)
  values_Q <- rep(1, nrow(K_QQ))
  K_QP <- outer(sample_ABC$x1, ref_ABC$x1, 
                Vectorize(function(x, y){return(exp(-sqrt(sum((x-y)^2))/0.1))}))
  score <- score_pdf(values_Q=values_Q, values_P=values_P,
                     K_PP=K_PP, K_QP=K_QP, K_QQ=K_QQ, 
                     discrete_Q =TRUE, discrete_P = FALSE)
  if(config$case=="geol"){
    #Compare against sample as well
    K_PP2 <- exp(-as.matrix(dist(samples_geol_ABC$x1, upper=TRUE, diag=TRUE))^2/0.1)
    K_QP2 <- outer(sample_ABC$x1, samples_geol_ABC$x1,
                   Vectorize(function(x, y){return(exp(-sqrt(sum((x-y)^2))/0.1))}))
    values_P2 <- rep(1, nrow(samples_geol_ABC))
    
    score2 <- score_pdf(values_Q=values_Q, values_P=values_P2,
                        K_PP=K_PP2, K_QP=K_QP2, K_QQ=K_QQ, 
                        discrete_Q =TRUE, discrete_P = TRUE)
    score_val2 <- c(score2)
  }
}
score_val <- c(score)

sizebatch <- 25
nbatch <- (nmax-nstart)/sizebatch
res <- data.frame(config, "PDF", t(c(score_val, rep(NA, nbatch))[1:(nbatch+1)]))
colnames(res)[-c(1:ncol(config))]<- c("reference", paste0("Score_", nstart+sizebatch*(0:nbatch)))
# 
# plot1 <- ggplot()+
#   geom_line(ref_ABC, mapping=aes(x=x1, y=ABCposterior))+
#   geom_histogram(data=sample_ABC, mapping=aes(x=x1, y=..density..))+
#   theme_bw()+
#   coord_cartesian(ylim=c(0, 7))
# plot2 <- ggplot()+
#   geom_hline(yintercept=10)+
#   geom_vline(xintercept=-1000)+
#   geom_line(x=seq(50, 1500, 25), 
#             y=score_val)+
#   theme_bw()+
#   coord_cartesian(xlim=c(50, 1500),
#                   ylim=c(0, 1))+
#   ylab("Score")+
#   xlab("n")
# ggarrange(plot1, plot2, nrow=2)
for(ibatch in seq(nbatch)){
  if(sampling=="RandomUnif"){
    loc_new <- runif(sizebatch)
    r_new <- rank(loc_new, ties.method = "first")
    loc_new <- sort(loc_new)
    source("../utils.R")
    tab_loc <- as.data.frame(table(loc_new))
    tab_loc$loc_new <- as.numeric(as.character(tab_loc$loc_new))
    if(config$case=="geol"){
      sample_new <- data.frame(t= unlist(sapply(seq(nrow(tab_loc)), function(i){
        sub_df <- df_ref[abs(round(tab_loc[i, 1], 14)- df_ref$x1)<=1/100, ]
        value_cdf <- interp2(x=unique(sub_df$x1),
                             y=u,
                             Z=matrix(sub_df$cdf, ncol=2),
                             yp=u, 
                             xp=rep(tab_loc[i, 1], length(u)))
        return(sample_from_custom_cdf(n=tab_loc[i, 2], t=u, cdf=value_cdf))
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
    loc_new <- runif(1)
    source("../utils.R")
    if(config$case=="geol"){
      sub_df <- df_ref[abs(round(loc_new, 14)- df_ref$x1)<=1/100, ]
      value_cdf <- interp2(x=unique(sub_df$x1),
                           y=u,
                           Z=matrix(sub_df$cdf, ncol=2),
                           yp=u, 
                           xp=rep(loc_new, length(u)))
      sample_new <- sample_from_custom_cdf(n=sizebatch, t=u, cdf=value_cdf)
      rm(sub_df)
    }else{
      source(paste0("../analytical_functions/ref_an_1D_f", config$fmed, "_noise", config$noise, ".R"))
      sample_new <-  sampling_values(unique_loc=loc_new, unique_freq=sizebatch)
    }
  }
  samples <- rbind(samples, sample_new)
  
  sample_ABC <- samples[samples$t <= treshold, ]
  if(nrow(sample_ABC)==0){
    score <- NA
    if(config$case=="geol"){
      score2 <- NA
      score_val2 <- c(score_val2, score2)
    }
  }else{
    K_QQ <- exp(-as.matrix(dist(sample_ABC$x1, upper=TRUE, diag=TRUE))^2/0.1)
    K_QP <- outer(sample_ABC$x1, ref_ABC$x1, 
                  Vectorize(function(x, y){return(exp(-sqrt(sum((x-y)^2))/0.1))}))
    values_Q <- rep(1, nrow(K_QQ))
    score <- score_pdf(values_Q=values_Q, values_P=values_P,
                       K_PP=K_PP, K_QP=K_QP, K_QQ=K_QQ, 
                       discrete_Q =TRUE, discrete_P = FALSE)
    if(config$case=="geol"){
      #Compare against sample as well
      K_QP2 <- outer(sample_ABC$x1, samples_geol_ABC$x1,
                     Vectorize(function(x, y){return(exp(-sqrt(sum((x-y)^2))/0.1))}))
      score2 <- score_pdf(values_Q=values_Q, values_P=values_P2,
                          K_PP=K_PP2, K_QP=K_QP2, K_QQ=K_QQ, 
                          discrete_Q =TRUE, discrete_P = TRUE)
      score_val2 <- c(score_val2, score2)
    }
  }
  score_val <- c(score_val, score)
  
  res <- data.frame(config, "PDF", t(c(score_val, rep(NA, nbatch))[1:(nbatch+1)]))
  colnames(res)[-c(1:ncol(config))]<- c("reference", paste0("Score_", nstart+sizebatch*(0:nbatch)))
  if(config$case=="geol"){
    temp <-data.frame(config, "sample", t(c(score_val2, rep(NA, nbatch))[1:(nbatch+1)]))
    colnames(temp)[-c(1:ncol(config))]<- c("reference", paste0("Score_", nstart+sizebatch*(0:nbatch)))
    res <-rbind(res,temp)
  }
  # 
  # plot1 <- ggplot()+
  #   geom_line(ref_ABC, mapping=aes(x=x1, y=ABCposterior))+
  #   geom_histogram(data=sample_ABC, mapping=aes(x=x1, y=..density..),
  #                  breaks=seq(0, 1,,31), alpha=0.1, col="black")+
  #   theme_bw()+
  #   coord_cartesian(ylim=c(0, 7))
  # xscore <- seq(50, 1500, 25)[1:length(score_val)]
  # indscore <- !is.na(score_val)
  # plot2 <- ggplot()+
  #   geom_hline(yintercept=10)+
  #   geom_vline(xintercept=-1000)+
  #   geom_line(data=data.frame(x=xscore[indscore], 
  #             y=score_val[indscore]),
  #             mapping=aes(x=x, y=y))+
  #   theme_bw()+
  #   coord_cartesian(xlim=c(50, 1500),
  #                   ylim=c(0, 1))+
  #   ylab("Score")+
  #   xlab("n")
  # show(ggarrange(plot1, plot2, nrow=2))
  write.csv(res,
            file=paste0("./res/scores/vanilla_sampling_", sampling, "_jobnumber_", job_number, ".txt"))
  
  write.csv(samples,
            file=paste0("./res/samples/vanilla_sampling_", sampling, "_jobnumber_", job_number, ".txt"))
  
}
