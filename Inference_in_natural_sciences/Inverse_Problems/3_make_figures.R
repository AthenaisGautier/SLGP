setwd("C:/Users/athen/Research/Codes/repo_PHD/Inference_in_natural_sciences/Inverse_Problems/")
# 
###Ref fields geol
library(ggplot2)
library(viridis)
library(ggpubr)
s <- 11
g <- 2
L <- list()
for(s in c(11, 25)){
  for(g in c(2, 162)){
    load(file=paste0("../ref_field_geol_", g, "_source_", s, ".RData"))
    treshold <- 0.15
    df_ref <- round(df_ref, 14)
    ref_ABC <- df_ref[df_ref$t==treshold, c("x1", "cdf")]
    ref_ABC$cdf <- ref_ABC$cdf/mean(ref_ABC$cdf)
    colnames(ref_ABC)[2] <- "ABCposterior"
    
    df_ref$type <- "Misfits and fitted field"
    samples_geol$type <- "Misfits and fitted field"
    ref_ABC$type <- "ABC posterior"
    L[[paste0(s, "-", g, "-1")]] <-  ggplot()+
      geom_raster(data=df_ref, mapping=aes(x=x1, y=t, fill=pdf), interpolate =TRUE)+
      geom_point(data=samples_geol, mapping=aes(x=x1, y=t), pch="x")+
      # facet_wrap(.~factor(type, levels=c("Misfits and fitted field", "ABC posterior")), ncol=1,
      #            scales = "free")+
      theme_bw()+
      scale_fill_gradientn(colors = c("white", viridis_pal(option="magma")(9)[10-1:8]), 
                           limits=c(0, 8), 
                           na.value = viridis_pal(option="magma")(9)[1],
                           guide=guide_colorbar(nrow=1, title="Probability density",
                                                barheight  = unit(2, units = "mm"),
                                                title.position = 'top',
                                                label.position = "bottom", 
                                                title.hjust = 0.5))+
      theme(legend.position="bottom", legend.direction = "horizontal",
            plot.title = element_text(hjust = 0.5))+
      xlab("Source depth (renormalized)")+
      ylab("Misfit value")+
      geom_hline(lty=2, col="blue", yintercept=0.15)+
      ggtitle(paste0("Geological structure ", g, "\nSource depth ~", round((s-1)/49*10, 2), "m. (",round((s-1)/49, 2), " renormalized)"))
    
    L[[paste0(s, "-", g, "-2")]] <-  ggplot()+
      geom_line(data=ref_ABC, mapping=aes(x=x1, y=ABCposterior))+
      theme_bw()+
      xlab("Source depth (renormalized)")+
      ylab("ABC posterior")+
      coord_cartesian(ylim=c(0, 8.25))
  }
}
ggarrange(L[["11-2-1"]], L[["25-2-1"]], 
          L[["11-2-2"]], L[["25-2-2"]], 
          L[["11-162-1"]], L[["25-162-1"]],
          L[["11-162-2"]], L[["25-162-2"]],
          common.legend = TRUE, legend = "bottom", ncol=2, nrow=4, heights = c(4, 2, 4, 2))
ggsave(file=paste0("./Figures/ref_geol_ABC.pdf"),
       width =8, height=6, scale=1.05)

###Scores
res <- read.csv(file="./res/combined_scores.txt")

temp <- colnames(res)
res<- res %>%
  pivot_longer(paste0("Score_", c(650, 675, 700)))%>%
  mutate(value=ifelse(geol==2 & source==11 & model=="SLGP" & strategy=="EIMAD", NA, value))%>%
  pivot_wider()
ind <- res$geol==2 & res$source==11 & res$model=="SLGP" & res$strategy=="EIMAD"
res <- res[, temp]

library(ggplot2)
library(viridis)
library(ggpubr)
library(tidyverse)

colnames(res)[8:66] <-unname( sapply(colnames(res)[8:66], function(x){
  strsplit(x, "Score_")[[1]][2]
}))
res <- res %>%
  dplyr::filter(!(reference %in% c("PDF VS mean", "sample VS mean")))%>%
  mutate(reference=ifelse(reference=="PDF VS prob", "PDF", reference))%>%
  mutate(reference=ifelse(reference=="sample VS prob", "sample", reference))

res2<- res %>%  
  pivot_longer(colnames(res)[8:66])%>%
  mutate(name=as.numeric(name))%>%
  rename(n=name)%>%
  group_by(case, source, geol, fmed, noise, reference, strategy, model, n)%>%
  summarise(value=median(value, na.rm = TRUE),
            finished_sim=n(),
            .groups="keep")%>%
  ungroup()%>%
  data.frame()

res2 %>%
  dplyr::filter(case=="geol")%>%
  dplyr::filter(reference=="PDF")%>%
  dplyr::filter(strategy!="RandomBatch")%>%
  dplyr::filter(model!="vanilla")%>%
  ggplot(mapping=aes(x=n, y=value, col=model, lty=strategy))+
  geom_line()+
  theme_bw()+
  facet_wrap(factor(paste0("Geological structure ", geol), levels=c(paste0("Geological structure ", c(2, 162))))~
               paste0("Source at ", round((source-1)/49, 2), " (normalised)"))+
  scale_y_log10()+
  ylab("Score value (in log scale)")+
  xlab("Sample size")
ggsave(file=paste0("./Figures/benchmark_geol.pdf"),
       width =6, height=5, scale=1.05)

res2 %>%
  dplyr::filter(case!="geol")%>%
  dplyr::filter(strategy!="RandomBatch")%>%
  dplyr::filter(model!="vanilla")%>%
  ggplot(mapping=aes(x=n, y=value, col=model, lty=strategy))+
  geom_line()+
  theme_bw()+
  facet_wrap(paste0("Median: f", fmed)~noise, ncol=4)+
  scale_y_log10()


config <- read.csv("config.txt")
config$jobnumber <- seq(nrow(config))
config %>%
  dplyr::filter(case=="geol" & geol==162 & source==25)


### Show some realizations !
g <- 162
s <- 25
job_number<-4

load(file=paste0("../ref_field_geol_", g, "_source_", s, ".RData"))
treshold <- 0.15
df_ref <- round(df_ref, 14)
ref_ABC <- df_ref[df_ref$t==treshold, c("x1", "cdf")]
ref_ABC$cdf <- ref_ABC$cdf/mean(ref_ABC$cdf)
colnames(ref_ABC)[2] <- "ABCposterior"

#Load a SLGP
sampling <- "RandomUnif"
res <- read.csv(paste0("./res/scores/SLGP_sampling_", sampling, "_jobnumber_", job_number, ".txt"))
res <- res[, -c(1)]
samples <- read.csv(file=paste0("./res/samples/SLGP_sampling_", sampling, "_jobnumber_", job_number, ".txt"))
samples <- samples[, -c(1)]
nrow(samples)

n <- 1500
load(file=paste0("./res/SLGPs/sampling_", sampling, "_jobnumber_", job_number, "_n_",
                 n, ".txt"))

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

par_basis_functions <- list(type="RFF",
                            nfreq=75,
                            dimension=2,
                            nu=5/2,
                            seed=config$rep)
initialize_basis_functions(par_basis_functions, verbose=1)
par_basis_functions$order_tot
name_index <-c("x1")

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

ggplot(df_plot%>%
                  dplyr::select(-starts_with("cdf"))%>%
                  pivot_longer(starts_with("pdf"))%>%
                  group_by(x1, t)%>%
                  summarise(value=mean(value), .groups="keep"), aes(x=x1, y=t))+
  geom_raster(mapping=aes(fill=value))+
  geom_point(data=samples)+
  geom_path(data=df_plot%>%
              dplyr::select(-starts_with("cdf"))%>%
              dplyr::filter(x1 %in% seq(0, 1, 0.1))%>%
              pivot_longer(starts_with("pdf"))%>%
              group_by(x1, t)%>%
              summarise(value=mean(value), .groups="keep")%>%
              group_by(x1)%>%
              mutate(value=value-min(value), .groups="keep")%>%
              data.frame()%>%
              arrange(t),
            mapping=aes(x=x1+value*0.05,
                y=t,
                group=x1), col="darkgrey")+
  geom_vline(xintercept = seq(0, 1, 0.1), lty=2, col="darkgrey")+
  theme_minimal()
mean_plot <-  df_ABC %>%
  pivot_longer(-x1)%>%
  group_by(x1)%>%
  summarise(value=mean(value))
df_ABC %>%
  pivot_longer(-x1)%>%
  ggplot()+
  geom_line(mapping=aes(x=x1, y=value, group=name), alpha=0.2)+
  geom_line(data=mean_plot,
            mapping=aes(x=x1, y=value),
            col="red")+
  coord_cartesian(ylim=c(0.2, 7))+
  geom_line(ref_ABC, mapping=aes(x=x1, y=ABCposterior), lty=2, col="blue")+
  ggtitle(n)+
  theme_bw()
