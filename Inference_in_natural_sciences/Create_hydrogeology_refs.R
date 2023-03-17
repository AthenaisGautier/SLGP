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



# To ensure the misfits are all in [0, 1] we get the highest simulated misfit and use it to normalise

if(TRUE){
  max_s <- 50
  max_g <- 2
  x <- paste(c("00000", max_g), collapse = "")
  x <- substring(x,nchar(x)-5+1)
  title_file <- paste(c("../datasets/Contaminant/distances/dist", x, "dt125z", max_s, ".txt"), collapse = "")
  samples <- data.frame(t=as.numeric(read.csv(title_file, header = FALSE)), x1 = c(t(matrix(rep(seq(0, 1, length.out = 50), 200), nrow=50))))
  max_val <- max(samples)
  rm(x, max_s, max_g, samples, title_file)
}


sources_considered <- c(11, 25)
geol_considered <- c(2, 162)

#Specify model
set.seed(1)
par_basis_functions <- list(type="RFF",
                            nfreq=250,
                            dimension=2,
                            nu=5/2,
                            seed=1)
initialize_basis_functions(par_basis_functions, verbose=1)
par_basis_functions$order_tot
name_index <-c("x1")

df_res <- data.frame(s=0, g=0, l_t=0, l_x=0, neglogpost=0)[c(0), ]
s <- 11
g <- 2

# for(s in sources_considered){
#   for(g in geol_considered){
#     x <- paste(c("00000", g), collapse = "")
#     x <- substring(x,nchar(x)-5+1)
#     title_file <- paste(c("../datasets/Contaminant/distances/dist", x, "dt125z", s, ".txt"), collapse = "")
#     samples <- data.frame(t=as.numeric(read.csv(title_file, header = FALSE)), 
#                           x1 = c(t(matrix(rep(seq(0, 1, length.out = 50), 200), nrow=50))))
#     samples$t <- samples$t/max_val
#     rm(x, title_file)
#     l_t <- l_x <- 0.05
#     for(l_t in c(seq(0.05, 0.15, 0.025), 0.2, .025)){
#       for(l_x in c(seq(0.05, 0.15, 0.025), 0.2, .025)){
#         starting_lengthscale <- c(l_t, l_x)
#         cat("G:", g, "S:", s, "l:", starting_lengthscale)
#         set.seed(1)
#         range_list <- heuristic_find_variance(par_basis_functions, 
#                                               lengthscale=starting_lengthscale,
#                                               nsimu=100, grid_size=51, plot=FALSE)
#         range_GP <- mean(range_list)
#         starting_sigma <- 5/range_GP
#         starting_point<- rnorm(par_basis_functions$order_tot)*starting_sigma
#         res_opt <- compute_MAP_without_hyperpar(samples, name_index,
#                                                 par_basis_functions,
#                                                 starting_point = starting_point,
#                                                 starting_sigma = starting_sigma,
#                                                 starting_lengthscale=starting_lengthscale,
#                                                 print_level=0)
#         df_plot <- evaluate_SLGP_grid(epsilon=res_opt$epsilon,
#                                       Xx=as.matrix(data.frame(x1=seq(0, 1,, 101))), 
#                                       name_index = name_index, 
#                                       par_basis_functions = par_basis_functions,
#                                       lengthscale = starting_lengthscale)
#         # ggplot(df_plot, aes(y=t, x=x1))+
#         #   geom_raster(mapping = aes(fill=pdf))+
#         #   # geom_point(data=samples)+
#         #   theme_bw()
#         df_res <- rbind(df_res, c(s, g, l_t, l_x, res_opt$value))
#         colnames(df_res) <- c("s", "g", "lt", "lx", "neglogpost")
#       }
#     }
#   }
# }
# write.csv(df_res, file="neg_log_post_for_geol_refs.txt", row.names=FALSE)
df_res <- read.csv(file="neg_log_post_for_geol_refs.txt")

df_res %>%
  dplyr::filter(lt <= 0.15 & lx  <= 0.15)%>%
  group_by(s, g)%>%
  mutate(neglogpost=100*neglogpost/min(neglogpost))%>%
  ggplot(aes(x=lx, y=lt, fill=neglogpost))+
  geom_tile()+
  facet_grid(paste0("Geological structure n°:", g)~paste0("Source at ", round(10*(s-1)/49, 1), "m"), scales = "free")+
  theme_bw()+
  geom_point(data=df_res%>%
               dplyr::group_by(s, g)%>%
               dplyr::filter(neglogpost==min(neglogpost))%>%
               mutate(neglogpost=100),
             pch="x", size=4)+
  scale_fill_viridis(option = "magma", direction = +1, 
                     name="Negative-log-posterior (in % of smallest value)",
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
  theme(legend.position="bottom", 
        legend.box = "horizontal")
ggsave("negative_log_posterior_profile.png", width=8, height=6)

starting_lengthscale <- df_res%>%
  dplyr::group_by(s, g)%>%
  dplyr::filter(neglogpost==min(neglogpost))%>%
  ungroup()%>%
  dplyr::select(lt, lx)
# All have the same optimum
starting_lengthscale <-as.numeric(starting_lengthscale[1, ])

# Redo Optim for reference fields.
for(s in sources_considered){
  for(g in geol_considered){
    x <- paste(c("00000", g), collapse = "")
    x <- substring(x,nchar(x)-5+1)
    title_file <- paste(c("../datasets/Contaminant/distances/dist", x, "dt125z", s, ".txt"), collapse = "")
    samples <- data.frame(t=as.numeric(read.csv(title_file, header = FALSE)), 
                          x1 = c(t(matrix(rep(seq(0, 1, length.out = 50), 200), nrow=50))))
    samples$t <- samples$t/max_val
    rm(x, title_file)
    starting_sigma <- 1
    set.seed(1)
    starting_point<- rnorm(par_basis_functions$order_tot)*starting_sigma
    res_opt <- compute_MAP_without_hyperpar(samples, name_index,
                                            par_basis_functions,
                                            starting_point = starting_point,
                                            starting_sigma = starting_sigma,
                                            starting_lengthscale=starting_lengthscale,
                                            print_level=0)
    df_plot <- evaluate_SLGP_grid(epsilon=res_opt$epsilon,
                                  Xx=as.matrix(data.frame(x1=seq(0, 1,, 101))), 
                                  name_index = name_index, 
                                  par_basis_functions = par_basis_functions,
                                  lengthscale = starting_lengthscale)
    ggplot(df_plot, aes(y=t, x=x1*10))+
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
      xlab("Source depth [m]")+
      ylab("Misfit value")+
      labs(caption=" \n ")
    ggsave(paste0("Figures/ref_field_geol_", g, "_source_", s, ".png"),
           width = 3.5, height=2, scale=1.75)
    ggplot(df_plot, aes(y=t, x=x1*10))+
      geom_raster(mapping = aes(fill=pdf),
                  interpolate=TRUE)+
      geom_point(data=samples, pch="x", col="grey20")+
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
      xlab("Source depth [m]")+
      ylab("Misfit value")+
      labs(caption=" \n ")
    ggsave(paste0("Figures/ref_field_geol_", g, "_source_", s, "_withdata.png"),
           width = 3.5, height=2, scale=1.75)
    treshold <- 0.15
    ggplot(df_plot, aes(y=t, x=x1*10))+
      geom_raster(mapping = aes(fill=pdf),
                  interpolate=TRUE)+
      geom_point(data=samples, pch="x", col="grey20")+
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
      xlab("Source depth [m]")+
      ylab("Misfit value")+
      geom_segment(x=0, xend=10, y=treshold, yend=treshold, lty=2, col="cornflowerblue")+
      labs(caption=paste0("Treshold at: ", treshold, "\n(", sum(samples$t<= treshold), " misfits below; ",
                          round(100*sum(samples$t<= treshold)/10000, 2), "% of simulations)"))
    ggsave(paste0("Figures/ref_field_geol_", g, "_source_", s, "_withdataandtreshold.png"),
           width = 3.5, height=2, scale=1.75)
    df_ABC <- df_plot %>%
      dplyr::filter(t==treshold)%>%
      mutate(x1=x1*10)%>%
      mutate(cdf=cdf/mean(cdf)/10)%>%
      dplyr::select(-pdf, -t)
    df_hist <- samples[samples$t <= treshold, ]
    df_hist$x1 <- df_hist$x1*10
    ggplot(df_ABC, aes(x=x1))+
      geom_line(aes(y=cdf))+
      theme_bw()+
      xlab("Source depth [m]")+
      ylab("ABC posterior")+
      geom_vline(xintercept=(s-1)/49*10, lty=2, col="midnightblue", lwd=1)+
      geom_histogram(df_hist, mapping=aes(y=..density..), 
                     breaks=seq(0-5/49, 10+5/49,,51), 
                     alpha=0.2, col="steelblue4", fill="steelblue")+
      labs(caption=paste0("Treshold at: ", treshold, "\n(", sum(samples$t<= treshold), " misfits below; ",
                          round(100*sum(samples$t<= treshold)/10000, 2), "% of simulations)"))
    ggsave(paste0("Figures/ref_field_geol_", g, "_source_", s, "_posteriorABC.png"),
           width = 3.5, height=2, scale=1.75)
    df_ABC <- df_plot %>%
      dplyr::filter(t==treshold)%>%
      mutate(x1=x1*10)%>%
      mutate(cdf=cdf/mean(cdf)/10)%>%
      dplyr::select(-pdf, -t)
    df_post <- df_plot %>%
      dplyr::filter(t==0)%>%
      mutate(x1=x1*10)%>%
      mutate(pdf=pdf/mean(pdf)/10)%>%
      dplyr::select(-cdf, -t)
    ggplot(df_post, aes(x=x1))+
      geom_line(aes(y=pdf))+
      theme_bw()+
      xlab("Source depth [m]")+
      ylab("ABC posterior")+
      geom_vline(xintercept=(s-1)/49*10, lty=2, col="midnightblue", lwd=1)+
      labs(caption=" \n ")
    ggsave(paste0("Figures/ref_field_geol_", g, "_source_", s, "_posterior.png"),
           width = 3.5, height=2, scale=1.75)
    df_ref <- df_plot
    par_basis_functions_ref <- par_basis_functions
    lengthscale_ref <- starting_lengthscale
    epsilon_ref <- res_opt$epsilon
    samples_geol <- samples
    save(par_basis_functions_ref, lengthscale_ref, epsilon_ref, df_ref, samples_geol,
         file=paste0("ref_field_geol_", g, "_source_", s, ".RData"))
  }
}
