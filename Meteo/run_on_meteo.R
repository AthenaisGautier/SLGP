setwd("C:/Users/athen/Research/Codes/repo_PHD/Meteo")

# Dataset
X <- read.csv("../datasets/data_meteo.txt", encoding = "latin1")
stations <- unique(X[, c(1, 8:13)])
stations$Station

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

library(grid)
library(gridExtra)
library(lemon)
library(ggrepel)
library(ggpubr)
library(metR)

packages_needed <- unique(packages_needed)
packages_needed <- packages_needed[!is.na(packages_needed)]
for(packname in packages_needed){
  library(packname, character.only=TRUE)
}
rm(file, to_source, packname, packages_needed, packs, txt)


world_map <- map_data("world", region="Switzerland")

range_temp <- c(-30, 40)
range_lat <- range(world_map$lat)
range_long <- range(world_map$long)
range_height <- c(0, 4810)
samples <- X[, c("tre200d0", "Latitude", "Longitude", "Station.height.m..a..sea.level")]
samples$tre200d0 <- (samples$tre200d0 - range_temp[1])/diff(range_temp)
samples$Latitude <- (samples$Latitude - range_lat[1])/diff(range_lat)
samples$Longitude <- (samples$Longitude - range_long[1])/diff(range_long)
samples$Station.height.m..a..sea.level <- (samples$Station.height.m..a..sea.level - range_height[1])/diff(range_height)
summary(samples)

colnames(samples) <- c("t", paste0("x", seq(3)))
summary(samples)
stations_rescaled <- unique(samples[, startsWith(colnames(samples), "x")])
name_index <- paste0("x", seq(3))


#Specify model
set.seed(1)
par_basis_functions <- list(type="RFF",
                            nfreq=250,
                            dimension=4,
                            nu=5/2,
                            seed=1)
# par_basis_functions <- list(type="regular",
#                             nfreq_x=2,
#                             nfreq_t=3,
#                             dimension=4,
#                             nu=5/2,
#                             seed=1)
initialize_basis_functions(par_basis_functions, verbose=1)
par_basis_functions$order_tot

starting_lengthscale <- 0.1
set.seed(1)
range_list <- heuristic_find_variance(par_basis_functions, lengthscale=starting_lengthscale,
                                      nsimu=100, grid_size=21, plot=TRUE)
range_GP <- mean(range_list)
starting_sigma <- 5/range_GP

# Show prior
# df_plot <- evaluate_SLGP_grid(epsilon=matrix(starting_sigma* rnorm(par_basis_functions$order_tot*10), nrow=10), 
#                               Xx=stations_rescaled, 
#                               name_index, 
#                               lengthscale=starting_lengthscale,
#                               par_basis_functions=par_basis_functions)
stations_rescaled$Station <- stations$Station
stations_rescaled$inBE <- stations$Canton=="BE"
# df_plot <- df_plot %>%
#   dplyr::select(-starts_with("cdf"))%>%
#   pivot_longer(-c("t", "x3"))

# df_plot <- merge(df_plot, stations_rescaled)
samples <- merge(samples, stations_rescaled)
samples <- samples[, c("t", name_index, "Station")]
samples$inBE <- samples$Station %in% stations$Station[stations$Canton=="BE"]
ggplot()+
  geom_histogram(data=samples, mapping=aes(x=t*diff(range_temp)+range_temp[1],
                                           y=..density..),
                 col="grey70", fill="grey95", size=0.5, 
                 breaks=seq(range_temp[1], range_temp[2], 2))+
  geom_density(data=samples, mapping=aes(x=t*diff(range_temp)+range_temp[1],
                                         y=..density..), lty=1,
               alpha=0, col="black")+
  # geom_line(mapping=aes(x=t*diff(range_temp)+range_temp[1], y=value/diff(range_temp), col=name),
  #           alpha=0.5)+
  facet_wrap(paste0(Station, "\n[", 
                    x3*diff(range_height)+range_height[1],
                    "m.]")~., ncol=5)+
  theme_bw()+
  ggtitle("Empirical distribution of temperatures")+
  theme(plot.title=element_text(hjust=0.5),
        legend.position = "none")+
  xlab("Temperature value [°C]")+
  ylab("Temperature distribution (marginalized over time)")+
  geom_rect(data = subset(stations_rescaled, stations_rescaled$inBE),
            inherit.aes = FALSE,
            fill = "cornflowerblue", col="blue", lty=2,
            alpha=0.2,
            xmin = -100, 
            xmax = 100,
            ymin = -1, 
            ymax = 1)
samples$Station <- NULL
stations_rescaled$Station <- NULL
stations_rescaled$inBE <- NULL

ggsave("./Figures/meteo_data.png", width=7, height=8)
ggsave("./Figures/meteo_data.pdf", width=7, height=8)

# ### Optim
# set.seed(1)
# starting_point<- rnorm(par_basis_functions$order_tot)*starting_sigma
# res_opt <- compute_MAP_without_hyperpar(samples, name_index,
#                                         par_basis_functions,
#                                         starting_point = starting_point,
#                                         starting_sigma = starting_sigma,
#                                         starting_lengthscale=starting_lengthscale,
#                                         print_level=0)
# 
# # Show MAP
# df_plot <- evaluate_SLGP_grid(epsilon=res_opt$epsilon, 
#                               Xx=stations_rescaled, 
#                               name_index, 
#                               lengthscale=starting_lengthscale,
#                               par_basis_functions=par_basis_functions)
# stations_rescaled$Station <- stations$Station
# 
# df_plot <- merge(df_plot, stations_rescaled)
# samples <- merge(samples, stations_rescaled)
# samples <- samples[, c("t", name_index, "Station")]
# 
# df_plot %>%
#   ggplot()+
#   geom_histogram(data=samples,
#                  mapping=aes(x=t*diff(range_temp)+range_temp[1], y=..density..), col="grey", alpha=0.1,
#                  breaks=seq(range_temp[1], range_temp[2],, 31))+
#   geom_line(mapping=aes(x=t*diff(range_temp)+range_temp[1], y=pdf/diff(range_temp)), col="red",
#             alpha=1)+
#   facet_wrap(Station~.)+
#   theme_bw()+
#   ggtitle("MAP")+
#   theme(legend.position = "none")+
#   xlab("Temperature value [°C]")+
#   ylab("Temperature distribution (marginalized over time)")
# df_plot %>%
#   ggplot(mapping=aes(x=t))+ 
#   stat_ecdf(data=samples, geom = "step")+
#   geom_line(mapping=aes(x=t, y=cdf), col="red",
#             alpha=1)+
#   facet_wrap(Station~.)+
#   theme_bw()+
#   ggtitle("MAP")+
#   theme(legend.position = "none")
# ggsave("MAP.png", width=8, height=8)
# samples$Station <- NULL
# stations_rescaled$Station <- NULL

#### 
# Now that it works: find the best lengthscale!
# 
list_lx <- seq(0.05, 0.6, 0.05)
list_lt <- seq(0.05, 0.3, 0.025)
# Lres <- list()
# df_res <- data.frame(lt=NA, lx1=NA, lx2=NA, lx3=NA, value=NA)[c(0), ]
# for(lx1 in list_lx){
#   lx2 <- lx1
#   for(lx3 in list_lx){
#     for(lt in list_lt){
#       set.seed(1)
#       starting_lengthscale <- c(lt, lx1, lx2, lx3)
#       names(starting_lengthscale) <- paste0("Lengthscale ", c("t", name_index))
#       print(starting_lengthscale)
#       range_list <- heuristic_find_variance(par_basis_functions, lengthscale=starting_lengthscale,
#                                             nsimu=100, grid_size=21, plot=FALSE)
#       range_GP <- mean(range_list)
#       starting_sigma <- 5/range_GP
#       # Optim
#       set.seed(1)
#       starting_point<- rnorm(par_basis_functions$order_tot)*starting_sigma
#       res_opt <- compute_MAP_without_hyperpar(samples, name_index,
#                                               par_basis_functions,
#                                               starting_point = starting_point,
#                                               starting_sigma = starting_sigma,
#                                               starting_lengthscale=starting_lengthscale,
#                                               print_level=0)
#       df_res <- rbind(df_res, c(starting_lengthscale, res_opt$value))
#       colnames(df_res) <- c(names(starting_lengthscale), "value")
#       if(res_opt$value==min(df_res$value)){
#         current_best <- res_opt
#       }
#       Lres[[length(Lres)+1]] <- res_opt
#       show(ggplot(df_res, aes(x=`Lengthscale t`, y=value))+
#         geom_line()+
#         facet_grid(`Lengthscale x1`~`Lengthscale x3`)+
#         theme_bw()+
#         geom_hline(yintercept=min(df_res$value), col="grey", lty=2))+
#         coord_cartesian(xlim=range(list_lt))
#     }
#   }
# }
# save.image("C:/Users/athen/Research/Codes/repo_PHD/run_on_meteo.RData")
# load("C:/Users/athen/Research/Codes/repo_PHD/run_on_meteo.RData")
# for(lx1 in seq(0.65, 0.8, 0.05)){
#   lx2 <- lx1
#   for(lx3 in seq(0.05, 0.35, 0.05)){
#     for(lt in seq(0.05, 0.2, 0.025)){
#       set.seed(1)
#       starting_lengthscale <- c(lt, lx1, lx2, lx3)
#       names(starting_lengthscale) <- paste0("Lengthscale ", c("t", name_index))
#       print(starting_lengthscale)
#       range_list <- heuristic_find_variance(par_basis_functions, lengthscale=starting_lengthscale,
#                                             nsimu=100, grid_size=21, plot=FALSE)
#       range_GP <- mean(range_list)
#       starting_sigma <- 5/range_GP
#       # Optim
#       set.seed(1)
#       starting_point<- rnorm(par_basis_functions$order_tot)*starting_sigma
#       res_opt <- compute_MAP_without_hyperpar(samples, name_index,
#                                               par_basis_functions,
#                                               starting_point = starting_point,
#                                               starting_sigma = starting_sigma,
#                                               starting_lengthscale=starting_lengthscale,
#                                               print_level=0)
#       df_res <- rbind(df_res, c(starting_lengthscale, res_opt$value))
#       colnames(df_res) <- c(names(starting_lengthscale), "value")
#       if(res_opt$value==min(df_res$value)){
#         current_best <- res_opt
#       }
#       Lres[[length(Lres)+1]] <- res_opt
#       show(ggplot(df_res, aes(x=`Lengthscale t`, y=value))+
#         geom_line()+
#         facet_grid(`Lengthscale x1`~`Lengthscale x3`)+
#         theme_bw()+
#         geom_hline(yintercept=min(df_res$value), col="grey", lty=2))+
#         coord_cartesian(xlim=range(list_lt))
#     }
#   }
# }
# save.image("C:/Users/athen/Research/Codes/repo_PHD/run_on_meteo.RData")
# # load("C:/Users/athen/Research/Codes/repo_PHD/run_on_meteo.RData")
# save(Lres, df_res, file="Meteo_df_res_and_Lres.RData")
load(file="Meteo_df_res_and_Lres.RData")

# Only display a subplot, for compacity purposes
ind <- df_res$`Lengthscale t` <= 0.2 &
  df_res$`Lengthscale x3`<= 0.3 &
  df_res$`Lengthscale x1`<= 0.5&
  df_res$`Lengthscale x1`>= 0.15
df_res <- df_res[ind, ]
Lres <- Lres[ind]
best_value <- range(df_res$value)
df_res <- df_res %>%
  mutate(`Lengthscale t`=`Lengthscale t`*100)%>%
  mutate(`Lengthscale x1`=ifelse(`Lengthscale x1`==0.05,
                                 paste0("Latitude\nlongitude\n05% of range"),
                                 paste0("Latitude\nlongitude\n", `Lengthscale x1`*100, "% of range")))%>%
  mutate(`Lengthscale x3`=ifelse(`Lengthscale x3`==0.05,
                                 paste0("Altitude\n05% of range"),
                                 paste0("Altitude\n", `Lengthscale x3`*100, "% of range")))%>%
  group_by_at(paste0("Lengthscale ", c(name_index)))%>%
  mutate(gives_min=any(value==best_value[1]))%>%
  mutate(to5percentmin=any(value<=0.99*best_value[1]))%>%
  ungroup()

### Plot for the optim profile
library(scales)
df_res %>%
  ggplot(aes(x=`Lengthscale t`, y=value))+
  geom_line()+
  geom_rect(data = subset(df_res, to5percentmin&!gives_min&`Lengthscale t`==10), 
            fill = "cornflowerblue", col="blue", lty=2,
            alpha=0.1,
            xmin = 100*(list_lt[1]-diff(range(list_lt))*10.04), 
            xmax = 100*(list_lt[length(list_lt)]+diff(range(list_lt))*10.04),
            ymin = 1*(best_value[1]-diff(best_value)*10.03), 
            ymax = 1*(best_value[2]+diff(best_value)*10.03))+
  geom_rect(data = subset(df_res, gives_min&`Lengthscale t`==10), 
            fill = "green", alpha=0.1,
            colour = "forestgreen", lty=1,
            xmin = 100*(list_lt[1]-diff(range(list_lt))*10.04), 
            xmax = 100*(list_lt[length(list_lt)]+diff(range(list_lt))*10.04),
            ymin = 1*(best_value[1]-diff(best_value)*10.03), 
            ymax = 1*(best_value[2]+diff(best_value)*10.03))+
  facet_grid(`Lengthscale x3`~`Lengthscale x1`)+
  theme_bw()+
  geom_hline(yintercept=min(df_res$value), col="grey")+
  xlab("Lengthscale for Temperature (in % of range)")+
  ylab("Negative log-posterior")+ 
  scale_y_continuous(labels = scientific)

ggsave("Figures/negative_log_likelihood_profile.png", width=8, height=6)
ggsave("Figures/negative_log_likelihood_profile.pdf", width=8, height=6)

samples$inBE <- NULL
# Optim without Bern

res_opt <- Lres[[which.min(df_res$value)]]
res_opt <- compute_MAP_without_hyperpar(samples[X$Canton != "BE", ], name_index,
                                        par_basis_functions,
                                        starting_point = res_opt$epsilon,
                                        starting_sigma = res_opt$sigma,
                                        starting_lengthscale=res_opt$lengthscale,
                                        print_level=0)

# Now: MCMC without Bern

Fit <- run_pCN(starting_point=res_opt$epsilon,
               samples[X$Canton != "BE", ], Iterations = 200000, Thinning = 2000, 
               par_basis_functions=par_basis_functions,
               name_index=name_index,
               starting_lengthscale = res_opt$lengthscale,
               starting_sigma = res_opt$sigma,
               n.approxInt = 101, 
               discretization.size = 101,
               beta=0.05, option="pCN")
Consort(Fit)

df_MCMC <- evaluate_SLGP_grid(epsilon=Fit$Posterior1, 
                              Xx=stations_rescaled, 
                              name_index, 
                              lengthscale=res_opt$lengthscale,
                              par_basis_functions=par_basis_functions)
df_MAP <- evaluate_SLGP_grid(epsilon=res_opt$epsilon, 
                             Xx=stations_rescaled, 
                             name_index, 
                             lengthscale=res_opt$lengthscale,
                             par_basis_functions=par_basis_functions)
stations_rescaled$Station <- stations$Station
stations_rescaled$inBE <- stations$Canton == "BE"
df_MAP <- merge(df_MAP, stations_rescaled)
df_MCMC <- merge(df_MCMC, stations_rescaled)
samples <- merge(samples, stations_rescaled)
samples <- samples[, c("t", name_index, "Station")]
samples$inBE <- ifelse(X$Canton == "BE", "In Bern (test)", "Not in Bern (train)")
df_MCMC$mean_pdf <- rowMeans(dplyr::select(df_MCMC, starts_with("pdf")))
df_MCMC$mean_cdf <- rowMeans(dplyr::select(df_MCMC, starts_with("cdf")))
df_plot <- df_MCMC[, c("t", "x1", "x2", "x3", "Station", "mean_pdf", "mean_cdf")]

df_MCMC %>%
  dplyr::select(c("t", "x1", "x2", "x3", "Station"), starts_with("pdf"))%>%
  pivot_longer(-c("t", "x1", "x2", "x3", "Station"))%>%
  ggplot()+
  geom_histogram(data=samples,
                 mapping=aes(x=t*diff(range_temp)+range_temp[1], y=..density..,
                             fill=inBE), col="grey", alpha=0.2,
                 breaks=seq(range_temp[1], range_temp[2],, 31))+
  geom_line(mapping=aes(x=t*diff(range_temp)+range_temp[1], y=value/diff(range_temp), group=name), col="cornflowerblue",
            alpha=0.1)+
  geom_line(data=df_plot, mapping=aes(x=t*diff(range_temp)+range_temp[1], y=mean_pdf/diff(range_temp)), 
            col="blue", lty=2)+
  geom_line(data=df_MAP, mapping=aes(x=t*diff(range_temp)+range_temp[1], y=pdf/diff(range_temp)), 
            col="red", lty=1)+
  # facet_wrap(Station~.)+
  theme_bw()+
  ggtitle("Posterior (MCMC)")+
  theme(legend.position = "bottom")+
  labs(fill="Localisation of the station")+
  xlab("Temperature value [°C]")+
  facet_wrap(paste0(Station, "\n[", 
                    x3*diff(range_height)+range_height[1],
                    "m.]")~., ncol=5)+
  scale_fill_manual(values=c( "cornflowerblue", "grey"))+
  scale_colour_manual(values=c("black", "blue"))+
  ylab("Temperature distribution (marginalized over time)")
# df_plot %>%
#   ggplot(mapping=aes(x=t))+ 
#   stat_ecdf(data=samples, geom = "step")+
#   geom_line(mapping=aes(x=t, y=cdf), col="red",
#             alpha=1)+
#   facet_wrap(Station~.)+
#   theme_bw()+
#   ggtitle("Posterior")+
#   theme(legend.position = "none")
ggsave("Posterior_full.png", width=8, height=8)
ggsave("Posterior_full.pdf", width=8, height=8)
samples$Station <- NULL
stations_rescaled$Station <- NULL
samples$inBE <- NULL
stations_rescaled$inBE <- NULL

mock_df <- data.frame(x=c(rep(-1000, 3), rep(-2000, 3)),
                      y=0, curve=rep(c("MAP estimate", "MCMC draws", "MCMC mean"), 2))

df_MCMC_rest <- df_MCMC%>%
  dplyr::select(starts_with("x") | t | starts_with("pdf") | Station)%>%
  pivot_longer(-c("t", "x1", "x2", "x3", "Station"))%>%
  mutate(t=t*diff(range_temp)+range_temp[1],
         value=value/diff(range_temp), 
         name=as.numeric(substr(name, 7, 11)))
ind_MCMC <- df_MCMC_rest$name %in% floor(seq(1, 100,, 10))

stations_rescaled$Station <- stations$Station
stations_rescaled$inBE <- stations$Canton == "BE"
ggplot(mock_df, mapping=aes(x=x, y=y, lty=curve, col=curve, alpha=curve))+
  geom_line()+
  theme_bw()+
  geom_histogram(data=samples,
                 mapping=aes(x=t*diff(range_temp)+range_temp[1], y=..density..), 
                 inherit.aes = FALSE,
                 col="grey70", fill="grey95", size=0.5, 
                 breaks=seq(range_temp[1], range_temp[2], 2))+
  facet_wrap(paste0(Station, "\n[", 
                    x3*diff(range_height)+range_height[1],
                    "m.]")~., ncol=5)+
  coord_cartesian(xlim=range_temp, ylim=c(0, 0.07))+
  scale_color_manual(values=c("red", "grey20", "blue"))+
  scale_linetype_manual(values=c(1, 1, 2))+
  scale_alpha_manual(values=c(1, 0.2, 1))+
  geom_line(data=data.frame(x=df_MCMC_rest$t,
                            y=df_MCMC_rest$value,
                            id=df_MCMC_rest$name, 
                            x3=df_MCMC_rest$x3,
                            Station=df_MCMC_rest$Station,
                            curve=mock_df$curve[2])[ind_MCMC, ],
            mapping=aes(group=id), alpha=0.5, lwd=0.3)+
  geom_line(data=data.frame(x=df_MAP$t*diff(range_temp)+range_temp[1],
                            y=df_MAP$pdf/diff(range_temp), 
                            x3=df_MAP$x3,
                            Station=df_MAP$Station,
                            curve=mock_df$curve[1]), lwd=0.5)+
  geom_line(data=data.frame(x=df_MCMC$t*diff(range_temp)+range_temp[1],
                            y=df_MCMC$mean_pdf/diff(range_temp),
                            curve=mock_df$curve[3],
                            x3=df_MCMC$x3,
                            Station=df_MCMC$Station), lwd=0.5)+
  xlab("Temperature value [°C]")+
  ylab("Temperature distribution\n(marginalized over time)")+
  theme(plot.title=element_text(size=10, hjust = 0.5),
        legend.position = "bottom",
        legend.title=element_blank())+
  geom_rect(data = subset(stations_rescaled, stations_rescaled$inBE),
            inherit.aes = FALSE,
            fill = "cornflowerblue", col="blue", lty=2,
            alpha=0.2,
            xmin = -100, 
            xmax = 100,
            ymin = -1, 
            ymax = 1)
ggsave("./Figures/posterior.png", width=7, height=8)
ggsave("./Figures/posterior.pdf", width=7, height=8)


library(readxl)
library(raster)
library(rgdal)
rel <- raster("../datasets/relief/DTM Switzerland, 50m, by Sonny/DTM Switzerland 50m.tif",
              crs="+proj=utm +zone=32 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs +type=crs")
# ,
# crs="+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333 +k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel +towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")
crs(rel)
rel <- projectRaster(rel, crs="+proj=longlat +datum=WGS84")
rel_spdf <- as(rel, "SpatialPixelsDataFrame")
rel <- as.data.frame(rel_spdf)


summary(rel)
colnames(rel) <- c("h", "long", "lat")
disc_long <- 101
disc_lat <- 101
rel2 <- rel
rel2$long <- round(rel2$long*(disc_long-1))/(disc_long-1)
rel2$lat <- round(rel2$lat*(disc_lat-1))/(disc_lat-1)
rel2 <- rel2 %>%
  group_by(long, lat) %>%
  summarise(h=median(h)) %>%
  data.frame()


world_map <- map_data("world", region="Switzerland")
rel2$long <- round(rel2$long*(disc_long-1))/(disc_long-1)
rel2$lat <- round(rel2$lat*(disc_lat-1))/(disc_lat-1)
# Plot beginning
switzerland <- ggplot() + 
  geom_raster(rel2, mapping=aes(x=long, y=lat, fill=h*0), fill="white") +
  geom_raster(rel2, mapping=aes(x=long, y=lat, fill=h), alpha=0.5) +
  theme_void()+
  theme(legend.position="bottom", 
        legend.box = "horizontal") +
  scale_fill_viridis(option = "turbo", direction = +1, 
                     name="Elevation [m]",
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
  xlab("Longitude") + ylab("Latitude") 
switzerland

Lin <- list()
samples$station <- X$Station
for(i in seq(nrow(stations))){
  Lin[[i]]<-list()
  Lin[[i]][["histogram"]] <- ggplot(X[X$Station==stations$Station[i], ])+
    theme_bw()+
    geom_histogram(mapping=aes(x=tre200d0,
                               y=..density..),
                   col="grey70", fill="grey95", size=0.5, 
                   breaks=seq(range_temp[1], range_temp[2], 2))+
    coord_cartesian(xlim=range_temp, ylim=c(0, 0.07))+
    xlab("Temperature value [°C]")+
    ylab("Temperature distribution\n(marginalized over time)")+
    ggtitle(paste0(stations$Station[i], " (", 
                   stations$Station.height.m..a..sea.level[i],
                   "m)")) + 
    theme(plot.title=element_text(size=10, hjust = 0.5))
  ind_MAP <- df_MAP$Station==stations$Station[i]
  ind_MCMC <- df_MCMC$Station==stations$Station[i]
  ind_MCMC2 <- df_MCMC_rest$Station==stations$Station[i] & df_MCMC_rest$name %in% seq(1, 100, 5)
  Lin[[i]][["data"]]<- ggplot(mock_df, mapping=aes(x=x, y=y, lty=curve, col=curve, alpha=curve))+
    geom_line()+
    theme_bw()+
    geom_histogram(data=X[X$Station==stations$Station[i], ],
                   inherit.aes=FALSE,
                   mapping=aes(x=tre200d0,
                               y=..density..), col="midnightblue",  fill="lightsteelblue", 
                   alpha=0.1,
                   breaks=seq(range_temp[1], range_temp[2], 2))+
    geom_line(data=data.frame(x=df_MCMC_rest$t[ind_MCMC2],
                              y=df_MCMC_rest$value[ind_MCMC2],
                              id=df_MCMC_rest$name[ind_MCMC2], 
                              curve=mock_df$curve[2]), mapping=aes(group=id), alpha=0.5)+
    geom_line(data=data.frame(x=df_MAP$t[ind_MAP]*diff(range_temp)+range_temp[1],
                              y=df_MAP$pdf[ind_MAP]/diff(range_temp),
                              curve=mock_df$curve[1]), lwd=1)+
    geom_line(data=data.frame(x=df_MCMC$t[ind_MCMC]*diff(range_temp)+range_temp[1],
                              y=df_MCMC$mean_pdf[ind_MCMC]/diff(range_temp),
                              curve=mock_df$curve[3]), lwd=1)+
    scale_color_manual(values=c("red", "darkgrey", "dodgerblue3"))+
    scale_linetype_manual(values=c(1, 1, 6))+
    scale_alpha_manual(values=c(1, 0.2, 1))+
    coord_cartesian(xlim=range_temp, ylim=c(0, 0.07))+
    xlab("Temperature value [°C]")+
    ylab("Temperature distribution\n(marginalized over time)")+
    ggtitle(paste0(stations$Station[i], " (", 
                   stations$Station.height.m..a..sea.level[i],
                   "m)")) + 
    theme(plot.title=element_text(size=10, hjust = 0.5),
          legend.position = "bottom",
          legend.title=element_blank())
}

ind_first_image <- c(3, 24, 13)
ind_stations <- c(3, 24, 13, 15, 27)
ind_in_be <- c(4, 12, 18)
plot_map <-switzerland+
  # theme_classic()+
  geom_point(data=stations[ind_first_image, ], mapping = aes(x=Longitude, y=Latitude)) +
  geom_label_repel(data=stations[ind_first_image, ],
                   aes(x = Longitude, y = Latitude, label=Station),
                   size=3, alpha=0.9, 
                   point.padding = 0., 
                   nudge_x = .25,
                   nudge_y = .25,
                   segment.curvature = -1e-20)+
  ggtitle("Stations location.") + 
  theme(plot.title=element_text(size=10, hjust = 0.5))
plot_map2 <-switzerland+
  # theme_classic()+
  geom_point(data=stations[ind_in_be, ], mapping = aes(x=Longitude, y=Latitude)) +
  geom_label_repel(data=stations[ind_in_be, ],
                   aes(x = Longitude, y = Latitude, label=Station),
                   size=3, alpha=0.9, 
                   point.padding = 0., 
                   nudge_x = .25,
                   nudge_y = .25,
                   segment.curvature = -1e-20)+
  ggtitle("Stations location.") + 
  theme(plot.title=element_text(size=10, hjust = 0.5))
plot_map3 <-switzerland+
  # theme_classic()+
  geom_point(data=stations[ind_stations, ], mapping = aes(x=Longitude, y=Latitude)) +
  geom_label_repel(data=stations[ind_stations, ],
                   aes(x = Longitude, y = Latitude, label=Station),
                   size=3, alpha=0.9, 
                   point.padding = 0., 
                   nudge_x = .25,
                   nudge_y = .25,
                   segment.curvature = -1e-20)+
  ggtitle("Stations location.") + 
  theme(plot.title=element_text(size=10, hjust = 0.5))

plot_grid1 <- grid.arrange(plot_map, 
                           Lin[[3]][["histogram"]], 
                           Lin[[13]][["histogram"]], 
                           Lin[[24]][["histogram"]],
                           ncol=2, nrow = 2)
ggsave(plot_grid1, file="./Figures/grid_with_data_noSLGP2.png", 
       width=8, height=2/3*8, scale=1)
ggsave(plot_grid1, file="./Figures/grid_with_data_noSLGP2.pdf", 
       width=8, height=2/3*8, scale=1)

legend1 <- get_legend(Lin[[1]][["data"]]+
                        theme(legend.title = element_text())+
                        labs(lty="Curve represented",
                             col="Curve represented",
                             alpha="Curve represented"))
legend2 <- get_legend(Lin[[1]][["data"]]+
                        theme(legend.position = "left",
                              legend.direction = "vertical",
                              legend.title = element_text())+
                        labs(lty="Curve represented",
                             col="Curve represented",
                             alpha="Curve represented"))
plot_grid2 <- grid.arrange(plot_map2, 
                           Lin[[4]][["data"]]+ 
                             theme(legend.position="none"), 
                           Lin[[12]][["data"]]+ 
                             theme(legend.position="none"), 
                           Lin[[18]][["data"]]+ 
                             theme(legend.position="none"),
                           legend2,
                           ncol=3, nrow = 2, 
                           widths = c(3, 3, 1.8), heights = c(3, 3),
                           layout_matrix = rbind(c(1, 2, 5),
                                                 c(3, 4, 5)))
ggsave(plot_grid2, file="./Figures/grid_without_data2.png", 
       width=9.5, height=5.5, scale=1/1.15)
ggsave(plot_grid2, file="./Figures/grid_without_data2.pdf", 
       width=9.5, height=5.5, scale=1/1.15)
plot_grid2 <- grid.arrange(plot_map2, 
                           Lin[[4]][["data"]]+ 
                             theme(legend.position="none"), 
                           Lin[[12]][["data"]]+ 
                             theme(legend.position="none"), 
                           Lin[[18]][["data"]]+ 
                             theme(legend.position="none"),
                           legend1,
                           ncol=3, nrow = 3, 
                           widths = c(3, 3, 3), heights = c(3, 3, 0.8),
                           layout_matrix = rbind(c(NA, 1, NA),
                                                 c(2, 3, 4),
                                                 c(5, 5, 5)))
ggsave(plot_grid2, file="./Figures/grid_without_data2_2.png", 
       width=9.5, height=5.5, scale=1/1.15)
ggsave(plot_grid2, file="./Figures/grid_without_data2_2.pdf", 
       width=9.5, height=5.5, scale=1/1.15)
plot_grid3 <- grid.arrange(plot_map3, 
                           Lin[[3]][["data"]]+ 
                             theme(legend.position="none"), 
                           Lin[[27]][["data"]]+ 
                             theme(legend.position="none"), 
                           Lin[[13]][["data"]]+ 
                             theme(legend.position="none"),
                           Lin[[24]][["data"]]+ 
                             theme(legend.position="none"),
                           Lin[[15]][["data"]]+ 
                             theme(legend.position="none"),
                           legend1,
                           ncol=3, nrow = 3, 
                           widths = c(3, 3, 3), heights = c(3, 3, 0.8),
                           layout_matrix = rbind(c(2, 1, 3),
                                                 c(4, 5, 6),
                                                 c(7, 7, 7)))
ggsave(plot_grid3, file="./Figures/grid_stbernard22.png", 
       width=9.5, height=5.5, scale=1/1.15)

ggsave(plot_grid3, file="./Figures/grid_stbernard22.pdf", 
       width=9.5, height=5.5, scale=1/1.15)


disc_long <- 101
disc_lat <- 101
df2 <- data.frame(rel2)
df2$lat <- (df2$lat - range_lat[1])/diff(range_lat)
df2$long <- (df2$long - range_long[1])/diff(range_long)
df2$h <- c((df2$h - range_height[1])/diff(range_height))
rownames(df2) <- NULL
summary(df2)
df2$long <- round(df2$long*(disc_long-1))/(disc_long-1)
df2$lat <- round(df2$lat*(disc_lat-1))/(disc_lat-1)
df2 <- df2 %>%
  group_by(long, lat) %>%
  summarise(h=median(h), .groups="keep") %>%
  data.frame()
colnames(df2) <- c("x2", "x1", "x3")
df2 <- df2[, paste0("x", seq(3))]

stations_discretized <- stations_rescaled
stations_discretized$x1 <- round(stations_discretized$x1*100)/100
stations_discretized$x2 <- round(stations_discretized$x2*100)/100
stations_discretized$x3 <- NULL
stations_discretized$inBE <- NULL

df_plot_meteo <- evaluate_SLGP_grid(epsilon=Fit$Posterior1, 
                                    Xx=df2, 
                                    name_index, 
                                    lengthscale=res_opt$lengthscale,
                                    par_basis_functions=par_basis_functions)

df_plot_meteo2 <- merge(df_plot_meteo, stations_discretized, by=c("x1", "x2"), all.x=TRUE)
df_plot_meteo2$subset2 <- abs(df_plot_meteo2$x1+1.125*df_plot_meteo2$x2-0.86)<=0.03
df_plot_meteo2$subset1 <- abs(df_plot_meteo2$x1+1.125*df_plot_meteo2$x2-0.86)<=0.005
plotmap <- switzerland +
  # geom_abline(slope=-1.125*diff(range_lat)/diff(range_long), 
  #             intercept = diff(range_lat)*
  #               (0.86+1.125*(range_long[1])/diff(range_long))+
  #               range_lat[1])+
  geom_segment(aes(x=0.175*diff(range_long)+range_long[1],
                   xend=0.68*diff(range_long)+range_long[1],
               y=47.12,
               yend=46.01), col="black", size=2)+
  geom_segment(aes(x=0.175*diff(range_long)+range_long[1],
                   xend=0.68*diff(range_long)+range_long[1],
                   y=47.12,
                   yend=46.01), col="grey", lty=1)
plotmap

df_plot_meteo2 <- subset(df_plot_meteo2, df_plot_meteo2$subset2)
# line y=-1.125 x+ 0.86
library(StereoMorph)
orth_vect <- orthogonalProjectionToLine(df_plot_meteo2[, c(1:2)], l1 = c(0, 0.86), 
                                        l2 = c(1, 0.86-1.125))
df_plot_meteo2$x1proj <- orth_vect[, 1]
df_plot_meteo2$x2proj <- orth_vect[, 2]
df_plot_meteo2$dproj <- sqrt((df_plot_meteo2$x1-df_plot_meteo2$x1proj)^2+
                               (df_plot_meteo2$x2-df_plot_meteo2$x2proj)^2)
df_plot_meteo2<- df_plot_meteo2%>%
  mutate(place_on_line = 1-(x1proj-min(x1proj))/diff(range(x1proj)))%>%
  group_by(dproj)

df_plot_meteo2_summary <- df_plot_meteo2 %>%
  dplyr::select(-x1proj, -x2proj, -Station)%>%
  rename(lat=x1,
         long=x2, 
         alt=x3,
         temp=t)%>%
  mutate(temp=temp*diff(range_temp)+range_temp[1],
         lat=lat*diff(range_lat)+range_lat[1],
         long=long*diff(range_long)+range_long[1],
         alt=alt*diff(range_height)+range_height[1])%>%
  pivot_longer(-c("lat", "long", "alt", "temp", "place_on_line", "dproj","subset1", "subset2"))%>%
  mutate(value=ifelse(substr(name, 1, 4)=="pdf_", value/diff(range_temp), value))%>%
  mutate(simu=substr(name, 5, 10),
         name=substr(name, 1, 3))%>%
  group_by(lat, long, alt, simu)%>%
  summarise(place_on_line=mean(place_on_line),
            subset1 = median(subset1),
            subset2=median(subset2),
            dproj=mean(dproj),
            mean = mean(temp*value*(name=="pdf"))*2,
            q20 = temp[which.min(1e10*(name=="pdf")+abs(value-0.2))],
            q30 = temp[which.min(1e10*(name=="pdf")+abs(value-0.3))],
            q10 = temp[which.min(1e10*(name=="pdf")+abs(value-0.1))],
            q90 = temp[which.min(1e10*(name=="pdf")+abs(value-0.9))],
            q40 = temp[which.min(1e10*(name=="pdf")+abs(value-0.4))],
            q60 = temp[which.min(1e10*(name=="pdf")+abs(value-0.6))],
            q80 = temp[which.min(1e10*(name=="pdf")+abs(value-0.8))],
            q50 = temp[which.min(1e10*(name=="pdf")+abs(value-0.5))],
            q70 = temp[which.min(1e10*(name=="pdf")+abs(value-0.7))], 
            .groups="keep")%>%
  ungroup()%>%
  data.frame()
stations_other <- subset(df_plot_meteo2,
                         !is.na(df_plot_meteo2$Station)&
                           df_plot_meteo2$t==0)
stations_other <- dplyr::select(stations_other, -starts_with("pdf"))
stations_other <- dplyr::select(stations_other, -starts_with("cdf"))
stations_other<- stations_other%>%
  rename(lat=x1,
         long=x2, 
         alt=x3,
         temp=t)%>%
  mutate(temp=temp*diff(range_temp)+range_temp[1],
         lat=lat*diff(range_lat)+range_lat[1],
         long=long*diff(range_long)+range_long[1],
         alt=alt*diff(range_height)+range_height[1])

plot1 <- 
  df_plot_meteo2_summary %>%
  filter(subset1==1)%>%
  dplyr::select(-mean)%>%
  pivot_longer(starts_with("q"))%>%
  mutate(name=paste0(substr(name, 2, 3), "%"))%>%
  rename(`Quantile level`=name)%>%
  group_by(lat, long, alt, place_on_line, `Quantile level`)%>%
  summarise(mean=mean(value),
            qL=quantile(value, 0.1),
            qH = quantile(value, 0.9))%>%
  ggplot(aes(x=place_on_line, y=alt))+
  geom_ribbon(mapping=aes(ymin=0, ymax=alt), alpha=0.2, col="black", lty=1)+
  geom_label_repel(data=stations_other, mapping=aes(label=Station))+
  theme_bw()+
  theme(legend.position="bottom", 
        legend.box = "horizontal")+
  ylab("Elevation above sea level [m]")+
  xlab("Position alongside the slice")
# geom_vline(data=data.frame(place_on_line=stations_other$place_on_line+
#                              c(-0.02, +0.01, +0.02, +0.01, 0, 0),
#                            alt=stations_other$alt),
#            mapping=aes(xintercept=place_on_line))
plot1

plot2 <- 
  df_plot_meteo2_summary %>%
  filter(subset1==1)%>%
  dplyr::select(-mean)%>%
  pivot_longer(starts_with("q"))%>%
  mutate(name=paste0(substr(name, 2, 3), "%"))%>%
  rename(`Quantile level`=name)%>%
  group_by(lat, long, alt, place_on_line, `Quantile level`)%>%
  summarise(mean=mean(value),
            qL=quantile(value, 0.1),
            qH = quantile(value, 0.9))%>%
  ggplot(aes(x=place_on_line, col=`Quantile level`, fill=`Quantile level`))+
  geom_line(mapping=aes(y=mean))+
  geom_ribbon(mapping=aes(ymin=qL, ymax=qH), alpha=0.2, lty=2)+
  theme_bw()+
  theme(legend.position="bottom", 
        legend.box = "horizontal")+
  ylab("Temperature value [°C]")+
  xlab("Position alongside the slice")+
  geom_vline(xintercept =stations_other$place_on_line+
               c(-0.02, +0.01, +0.02, +0.01, 0, 0), lty=2)+
  geom_label(data=stations_other, mapping=aes(label=Station, 
                                              y=c(-22, -17, -22,-12,-17,-22), x=place_on_line), inherit.aes = FALSE)
plot2 


plot_grid <- grid.arrange(plotmap+ theme(aspect.ratio=0.5), 
                           plot1, plot2,
                           ncol=2, nrow = 2, 
                           widths = c(3, 3), heights = c(3, 4),
                           layout_matrix = rbind(c(1, 1),
                                                 c(2, 3)),
                          respect=FALSE)
ggsave(plot_grid, file="./Figures/slice_switzerland.png", 
       width=9.5, height=5.5, scale=1.1)
ggsave(plot_grid, file="./Figures/slice_switzerland.pdf", 
       width=9.5, height=5.5, scale=1.1)

subset <- df_plot_meteo2$x1
df_plot_meteo_summary <- df_plot_meteo %>%
  rename(lat=x1,
         long=x2, 
         alt=x3,
         temp=t)%>%
  mutate(temp=temp*diff(range_temp)+range_temp[1],
         lat=lat*diff(range_lat)+range_lat[1],
         long=long*diff(range_long)+range_long[1],
         alt=alt*diff(range_height)+range_height[1])%>%
  pivot_longer(-c("lat", "long", "alt", "temp"))%>%
  mutate(value=ifelse(substr(name, 1, 4)=="pdf_", value/diff(range_temp), value))%>%
  mutate(simu=substr(name, 5, 10),
         name=substr(name, 1, 3))%>%
  group_by(lat, long, alt, simu)%>%
  summarise(mean = mean(temp*value*(name=="pdf"))*2,
            q50 = temp[which.min(1e10*(name=="pdf")+abs(value-0.5))],
            q10 = temp[which.min(1e10*(name=="pdf")+abs(value-0.1))],
            q90 = temp[which.min(1e10*(name=="pdf")+abs(value-0.9))],
            q25 = temp[which.min(1e10*(name=="pdf")+abs(value-0.25))],
            q75 = temp[which.min(1e10*(name=="pdf")+abs(value-0.75))],
            inv25 = value[which.min(abs(temp-25)+1e10*(name=="pdf"))],
            inv20 = value[which.min(abs(temp-20)+1e10*(name=="pdf"))],
            inv15 = value[which.min(abs(temp-15)+1e10*(name=="pdf"))],
            inv10 = value[which.min(abs(temp-10)+1e10*(name=="pdf"))],
            invminus10 = value[which.min(abs(temp+10)+1e10*(name=="pdf"))],
            inv5 = value[which.min(abs(temp-5)+1e10*(name=="pdf"))],
            invminus5 = value[which.min(abs(temp+5)+1e10*(name=="pdf"))],
            inv = value[which.min(abs(temp+0)+1e10*(name=="pdf"))], .groups="keep")%>%
  ungroup()%>%
  data.frame()

df_plot_meteo_summary <- df_plot_meteo_summary %>%
  ungroup()%>%
  pivot_longer(-c("lat", "long", "alt", "simu"))%>%
  group_by(lat, long, alt, name)%>%
  summarise(expectation=mean(value),
            sd=sd(value), .groups="keep")%>%
  ungroup()%>%
  rename(quantity=name)%>%
  pivot_longer(-c("lat", "long", "alt", "quantity"))%>%
  
  pivot_wider(names_from = quantity, values_from = value)%>%
  ungroup()%>%
  data.frame()
ind_exp <- df_plot_meteo_summary$name=="expectation"
plot1 <- ggplot(df_plot_meteo_summary[ind_exp, ],
                aes(x=long, y=lat, fill=q50))+
  geom_raster()+
  theme_void()+
  theme(legend.position="bottom", 
        legend.box = "horizontal") +
  scale_fill_viridis(option = "magma", direction = -1, 
                     name="Expected value [°C]",
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
  xlab("Longitude") + ylab("Latitude") 
plot2 <- ggplot(df_plot_meteo_summary[!ind_exp, ],
                aes(x=long, y=lat, fill=q50))+
  geom_raster()+
  theme_void()+
  theme(legend.position="bottom", 
        legend.box = "horizontal") +
  scale_fill_viridis(option = "magma", direction = -1, 
                     name="Standard deviation [°C]",
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
  xlab("Longitude") + ylab("Latitude") 
plot_grid <- grid.arrange(plot1, plot2, ncol=2)
ggsave(plot_grid, file="./Figures/grid_about_q50.png", 
       width=9, height=3.5, scale=1)
ggsave(plot_grid, file="./Figures/grid_about_q50.pdf", 
       width=9, height=3.5, scale=1)

plot1 <- ggplot(df_plot_meteo_summary[ind_exp, ],
                aes(x=long, y=lat, fill=mean))+
  geom_raster()+
  theme_void()+
  theme(legend.position="bottom", 
        legend.box = "horizontal") +
  scale_fill_viridis(option = "magma", direction = -1, 
                     name="Expected value [°C]",
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
  xlab("Longitude") + ylab("Latitude") 
plot2 <- ggplot(df_plot_meteo_summary[!ind_exp, ],
                aes(x=long, y=lat, fill=mean))+
  geom_raster()+
  theme_void()+
  theme(legend.position="bottom", 
        legend.box = "horizontal") +
  scale_fill_viridis(option = "magma", direction = -1, 
                     name="Standard deviation [°C]",
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
  xlab("Longitude") + ylab("Latitude") 
plot_grid <- grid.arrange(plot1, plot2, ncol=2)
ggsave(plot_grid, file="./Figures/grid_about_mean.png", 
       width=9, height=3.5, scale=1)
ggsave(plot_grid, file="./Figures/grid_about_mean.pdf", 
       width=9, height=3.5, scale=1)

