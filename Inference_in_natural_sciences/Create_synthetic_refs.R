setwd("C:/Users/athen/Research/Codes/repo_PHD/Inference_in_natural_sciences")

# Functions and packages
to_source <- list.files(path = "../functions_to_source/")
packages_needed <- c()
for(file in to_source){
  source(paste0("../functions_to_source/", file))
  library(readr)
  txt <- read_file(paste0("../functions_to_source/", file))
  library(qdapRegex)
  packs <- rm_between(txt, c("library(", "require("), c(")", ")"), extract=TRUE)
  packages_needed<- c(packages_needed, packs[[1]])
}
library(conflicted)
conflict_prefer("dinvchisq", "LaplacesDemon") 
conflict_prefer("dinvgamma", "LaplacesDemon") 
conflict_prefer("rinvchisq", "LaplacesDemon") 
conflict_prefer("rinvgamma", "LaplacesDemon") 

packages_needed <- unique(packages_needed)
packages_needed <- packages_needed[!is.na(packages_needed)]
for(packname in packages_needed){
  library(packname, character.only=TRUE)
}
rm(file, to_source, packname, packages_needed, packs, txt)

starting_lengthscale <- c(0.1, 0.1)
starting_sigma <- 1
# Do Optim for reference fields.
for(med in c(1, 2)){
  for(noise in seq(4)){
    source(paste0("./analytical_functions/ref_an_1D_f", med, "_noise", noise, ".R"))
    X <- expand.grid(seq(0, 1,, 101), seq(0, 1,, 101))
    colnames(X) <- c("t", "x1")
    df_plot <- ref_field(X)
    ggplot(df_plot, aes(y=t, x=x1))+
      geom_raster(mapping = aes(fill=pdf),
                  interpolate=TRUE)+
      # geom_point(data=samples)+
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
      xlab("Index x")+
      ylab("Response value")+
      labs(caption=" \n ")
    ggsave(paste0("Figures/ref_field_f", med, "_noise_", noise, ".png"),
           width = 3.5, height=2, scale=1.75)
    treshold <- 0.15
    ggplot(df_plot, aes(y=t, x=x1))+
      geom_raster(mapping = aes(fill=pdf),
                  interpolate=TRUE)+
      # geom_point(data=samples)+
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
      xlab("Index x")+
      ylab("Response value")+
      geom_segment(x=0, xend=10, y=treshold, yend=treshold, lty=2, col="cornflowerblue")+
      labs(caption=paste0("Treshold at: ", treshold, "\n"))
    ggsave(paste0("Figures/ref_field_f", med, "_noise_", noise,"_withtreshold.png"),
           width = 3.5, height=2, scale=1.75)
    df_ABC <- df_plot %>%
      dplyr::filter(t==treshold)%>%
      mutate(x1=x1)%>%
      mutate(cdf=cdf/mean(cdf))%>%
      dplyr::select(-pdf, -t)
    ggplot(df_ABC, aes(x=x1))+
      geom_line(aes(y=cdf))+
      theme_bw()+
      xlab("Index x")+
      ylab("Response value")+
      labs(caption=paste0("Treshold at: ", treshold, "\n"))
    ggsave(paste0("Figures/ref_field_f", med, "_noise_", noise, "_posteriorABC.png"),
           width = 3.5, height=2, scale=1.75)
    
  }
}
