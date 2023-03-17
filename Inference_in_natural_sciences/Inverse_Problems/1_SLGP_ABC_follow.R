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

print(paste0("SLGP, sampling ", sampling, " job ", job_number))

# Functions and packages
to_source <- list.files(path = "../../functions_to_source/")
packages_needed <- c()
for(file in to_source){
  source(paste0("../../functions_to_source/", file))
  library(readr)
  txt <- read_file(paste0("../../functions_to_source/", file))
  library(qdapRegex)
  packs <- rm_between(txt, c("library(", "require("), c(")", ")"), extract=TRUE)
  packages_needed<- c(packages_needed, packs[[1]])
}
library(conflicted)
conflict_prefer("dinvchisq", "LaplacesDemon") 
conflict_prefer("dinvgamma", "LaplacesDemon") 
conflict_prefer("rinvchisq", "LaplacesDemon") 
conflict_prefer("rinvgamma", "LaplacesDemon") 
library(Hmisc)
conflict_prefer("summarize", "dplyr") 

packages_needed <- unique(packages_needed)
packages_needed <- packages_needed[!is.na(packages_needed)]
for(packname in packages_needed){
  library(packname, character.only=TRUE)
}
rm(file, to_source, packname, packages_needed, packs, txt)


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


par_basis_functions <- list(type="RFF",
                            nfreq=75,
                            dimension=2,
                            nu=5/2,
                            seed=config$rep)
initialize_basis_functions(par_basis_functions, verbose=1)
par_basis_functions$order_tot
name_index <-c("x1")
starting_lengthscale <- c(0.15, 0.15)
range_list <- heuristic_find_variance(par_basis_functions, lengthscale=starting_lengthscale,
                                      nsimu=100, grid_size=51, plot=FALSE)
range_GP <- mean(range_list)
starting_sigma <- 5/range_GP

res <- read.csv(paste0("./res/scores/SLGP_sampling_", sampling, "_jobnumber_", job_number, ".txt"))
res <- res[, -c(1)]
samples <- read.csv(file=paste0("./res/samples/SLGP_sampling_", sampling, "_jobnumber_", job_number, ".txt"))
samples <- samples[, -c(1)]

load(file=paste0("./res/SLGPs/sampling_", sampling, "_jobnumber_", job_number, "_n_",
                 nrow(samples), ".txt"))

sizebatch <- 25
nbatch <- (nmax-nstart)/sizebatch
nsimu <- 100
batchstart <- (nrow(samples)-nstart)/sizebatch+1
ibatch <- batchstart
if(nrow(samples) != nmax){
  
  for(ibatch in seq(batchstart, nbatch)){
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
      mean_cdf <- rowMeans(dplyr::select(df_plot, starts_with("cdf")))
      i <- 1
      crit <- matrix(NA, ncol=length(candidates), nrow=nsimu)
      for(i in seq_along(candidates)){
        print(paste0("Simu for candidate ", i, "/", length(candidates)))
        xcand <- candidates[i]
        ind <- abs(df_plot$x1-xcand) <= 1e-10
        groups <- c(sapply(seq(nsimu), function(j){rep(j, sizebatch)}))
        fake_samp <- sample_from_custom_cdf(n=nsimu*sizebatch, t=df_plot$t[ind],
                                            cdf=mean_cdf[ind])
        fake_samp_likes <- sapply(seq(100), function(j){
          return(approx(x=df_plot$t[ind], y=df_plot[ind, (j-1)*2+3],
                        xout=fake_samp)$y)
        })
        weights <- apply(fake_samp_likes, 2, 
                         FUN=function(c){
                           tapply(c, INDEX=groups, FUN=function(x){
                             return(exp(sum(log(x))))
                           })
                         })
        weights <- weights/rowSums(weights)
        IV <- sapply(seq(nsimu), function(j){
          mean(sapply(seq(nrow(df_ABC)), function(k){
            wvars <- wtd.var(x=df_ABC[k, -c(1)], weights = weights[j, ], normwt = TRUE)
          }))
        })
        crit[, i] <- IV
      }
      
      # plot(candidates[c(sapply(seq(length(candidates)), function(i){rep(i, nsimu)}))], c(crit), pch="x")
      model <- loess(crit~x,
                     data=data.frame(x=candidates[c(sapply(seq(length(candidates)), function(i){rep(i, nsimu)}))],
                                     crit=c(-crit)), span = 0.15)
      # lines(seq(0, 1,, 101), predict(model, seq(0, 1,, 101)))
      loc_new <- candidates[which.min(predict(model, candidates))]
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
      mean_cdf <- rowMeans(dplyr::select(df_plot, starts_with("cdf")))
      i <- 1
      crit <- matrix(NA, ncol=length(candidates), nrow=nsimu)
      for(i in seq_along(candidates)){
        xcand <- candidates[i]
        ind <- abs(df_plot$x1-xcand) <= 1e-10
        groups <- c(sapply(seq(nsimu), function(j){rep(j, sizebatch)}))
        fake_samp <- sample_from_custom_cdf(n=nsimu*sizebatch, t=df_plot$t[ind],
                                            cdf=mean_cdf[ind])
        fake_samp_likes <- sapply(seq(100), function(j){
          return(approx(x=df_plot$t[ind], y=df_plot[ind, (j-1)*2+3],
                        xout=fake_samp)$y)
        })
        weights <- apply(fake_samp_likes, 2, 
                         FUN=function(c){
                           tapply(c, INDEX=groups, FUN=function(x){
                             return(exp(sum(log(x))))
                           })
                         })
        weights <- weights/rowSums(weights)
        
        IMAD <- sapply(seq(nsimu), function(j){
          mean(sapply(seq(nrow(df_ABC)), function(k){
            wtdmed <- wtd.quantile(x=as.numeric(unname(df_ABC[k, -c(1)])), 
                                   weights = weights[j, ], 
                                   normwt = TRUE,
                                   type='quantile',
                                   probs=c(0.5))
            MAD <- wtd.mean(x=abs(df_ABC[k, -c(1)]-wtdmed), weights = weights[j, ])
            return(MAD)
          }))
        })
        crit[, i] <- IMAD
      }
      # plot(candidates[c(sapply(seq(length(candidates)), function(i){rep(i, nsimu)}))], c(crit), pch="x")
      model <- loess(crit~x,
                     data=data.frame(x=candidates[c(sapply(seq(length(candidates)), function(i){rep(i, nsimu)}))],
                                     crit=c(-crit)), span = 0.15)
      # lines(seq(0, 1,, 101), predict(model, seq(0, 1,, 101)))
      loc_new <- candidates[which.min(predict(model, candidates))]
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
    
    samples <- rbind(samples, sample_new)
    
    
    starting_point<- Fit$Posterior1[nrow(Fit$Posterior1), ]
    # beta <- 0.905*sqrt(50/nrow(samples))-0.11
    Fit <- run_pCN(starting_point=starting_point,
                   samples, Iterations =10000, Thinning = 1, 
                   par_basis_functions=par_basis_functions,
                   name_index=name_index,
                   starting_lengthscale = starting_lengthscale,
                   starting_sigma = starting_sigma,
                   n.approxInt = 101, 
                   discretization.size = 101,
                   beta=beta, option="pCN")
    # plot(Fit$Posterior1[, 1])
    # plot(Fit$Posterior1[5000+seq(1, 5000, 50), 1])
    Acceptance.Rate <- AcceptanceRate(Fit$Posterior1[-seq(1, 5000), 1:2])[1]
    if(abs(Acceptance.Rate-0.25)>0.05){
      while(abs(Acceptance.Rate-0.25)>0.05){
        cat(beta, Acceptance.Rate, "\n")
        if(Acceptance.Rate>0.25){
          beta <- min(1.025*beta, 1)
        }else{
          beta <- 0.975*beta
        }
        starting_point<- Fit$Posterior1[nrow(Fit$Posterior1), ]
        Fit <- run_pCN(starting_point=starting_point,
                       samples, Iterations =10000, Thinning = 1, 
                       par_basis_functions=par_basis_functions,
                       name_index=name_index,
                       starting_lengthscale = starting_lengthscale,
                       starting_sigma = starting_sigma,
                       n.approxInt = 101, 
                       discretization.size = 101,
                       beta=beta, option="pCN")
        Acceptance.Rate <- AcceptanceRate(Fit$Posterior1[-seq(1, 5000), 1:2])[1]
      }
    }
    
    df_plot <- evaluate_SLGP_grid(epsilon=Fit$Posterior1[5000+seq(1, 5000, 50), ], 
                                  Xx=seq(0, 1,, 101), 
                                  name_index, 
                                  lengthscale = starting_lengthscale,
                                  par_basis_functions=par_basis_functions)
    
    df_ABC <- df_plot %>%
      dplyr::filter(t==treshold)%>%
      dplyr::select(-starts_with("pdf"), -t)%>%
      pivot_longer(-x1)%>%
      group_by(name)%>%
      mutate(value=value/mean(value))%>%
      pivot_wider()%>%
      data.frame()
    
    # 
    # plot1 <- ggplot(df_plot%>%
    #                   dplyr::select(-starts_with("cdf"))%>%
    #                   pivot_longer(starts_with("pdf"))%>%
    #                   group_by(x1, t)%>%
    #                   summarise(value=mean(value), .groups="keep"), aes(x=x1, y=t))+
    #   geom_raster(mapping=aes(fill=value))+
    #   geom_point(data=samples)+
    #   geom_path(data=df_plot%>%
    #               dplyr::select(-starts_with("cdf"))%>%
    #               dplyr::filter(x1 %in% seq(0, 1, 0.1))%>%
    #               pivot_longer(starts_with("pdf"))%>%
    #               group_by(x1, t)%>%
    #               summarise(value=mean(value), .groups="keep")%>%
    #               group_by(x1)%>%
    #               mutate(value=value-min(value), .groups="keep")%>%
    #               data.frame()%>%
    #               arrange(t),
    #             mapping=aes(x=x1+value*0.05,
    #                         y=t,
    #                         group=x1), col="darkgrey")+
    #   geom_vline(xintercept = seq(0, 1, 0.1), lty=2, col="darkgrey")+
    #   theme_minimal()
    # # 
    # mean_plot <-  df_ABC %>%
    #   pivot_longer(-x1)%>%
    #   group_by(x1)%>%
    #   summarise(value=mean(value))
    # plot2 <- df_ABC %>%
    #   pivot_longer(-x1)%>%
    #   ggplot()+
    #   geom_line(mapping=aes(x=x1, y=value, group=name))+
    #   geom_line(data=mean_plot, 
    #             mapping=aes(x=x1, y=value),
    #             col="red")+
    #   coord_cartesian(ylim=c(0.2, 1.2))+
    #   ggtitle(nrow(samples))+
    #   theme_bw()
    # show(grid.arrange(plot1, plot2, nrow=2))
    
    # The probabilistic pred
    values_Q <- df_ABC[, -c(1)]
    if(ibatch==batchstart){
      K_QQ <- K_QP <- K_PP <- exp(-as.matrix(dist(df_ABC$x1, upper=TRUE, diag=TRUE))^2/0.1)
      values_P <- ref_ABC$ABCposterior
      
      score <- score_pdf(values_Q=values_Q, values_P=values_P, 
                         K_PP=K_PP, K_QP=K_QP, K_QQ=K_QQ, 
                         discrete_Q =FALSE, discrete_P = FALSE)
      
      score_val <- c(unname(t(res[1, 7+1:batchstart])), score)
      
      if(config$case=="geol"){
        K_PP2 <- exp(-as.matrix(dist(samples_geol_ABC$x1, upper=TRUE, diag=TRUE))^2/0.1)
        K_QP2 <- outer(df_ABC$x1, samples_geol_ABC$x1,
                       Vectorize(function(x, y){return(exp(-sqrt(sum((x-y)^2))/0.1))}))
        values_P2 <- rep(1, nrow(samples_geol_ABC))
        score2 <- score_pdf(values_Q=values_Q, values_P=values_P2, 
                            K_PP=K_PP2, K_QP=K_QP2, K_QQ=K_QQ, 
                            discrete_Q =FALSE, discrete_P = TRUE)
        score_val2 <- c(unname(t(res[2, 7+1:batchstart])), score2)
      }
    }else{
      score <- score_pdf(values_Q=values_Q, values_P=values_P, 
                         K_PP=K_PP, K_QP=K_QP, K_QQ=K_QQ, 
                         discrete_Q =FALSE, discrete_P = FALSE)
      
      # plot(df_ABC, type="l")
      score_val <- c(score_val, score)
      if(config$case=="geol"){
        score2 <- score_pdf(values_Q=values_Q, values_P=values_P2, 
                            K_PP=K_PP2, K_QP=K_QP2, K_QQ=K_QQ, 
                            discrete_Q =FALSE, discrete_P = TRUE)
        score_val2 <- c(score_val2, score2)
      }
    }
    
    
    
    res <- data.frame(config, "PDF", t(c(score_val, rep(NA, nbatch))[1:(nbatch+1)]))
    colnames(res)[-c(1:ncol(config))]<- c("reference", paste0("Score_", nstart+sizebatch*(0:nbatch)))
    if(config$case=="geol"){
      temp <-data.frame(config, "sample", t(c(score_val2, rep(NA, nbatch))[1:(nbatch+1)]))
      colnames(temp)[-c(1:ncol(config))]<- c("reference", paste0("Score_", nstart+sizebatch*(0:nbatch)))
      res <-rbind(res,temp)
    }
    
    write.csv(res,
              file=paste0("./res/scores/SLGP_sampling_", sampling, "_jobnumber_", job_number, ".txt"))
    save(beta, Fit,
         file=paste0("./res/SLGPs/sampling_", sampling, "_jobnumber_", job_number, "_n_",
                     nrow(samples), ".txt"))
    write.csv(samples,
              file=paste0("./res/samples/SLGP_sampling_", sampling, "_jobnumber_", job_number, ".txt"))
  }
  
  
}