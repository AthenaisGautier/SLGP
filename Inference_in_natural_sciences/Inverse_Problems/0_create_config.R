setwd("C:/Users/athen/Research/Codes/repo_PHD/Inference_in_natural_sciences/Inverse_Problems")

test_cases <- rbind(expand.grid(c("geol"), c(11, 25), c(2, 162), NA, NA),
                    expand.grid(c("analytical"), NA, NA, c(1, 2), c(1, 2, 3, 4)))
rep <- seq(50)

config <- expand.grid(seq(nrow(test_cases)), rep)
config <- cbind(test_cases[config[, 1], ], config[, 2])
config$Var1 <- as.character(config$Var1)
colnames(config) <- c("case", "source", "geol", "fmed", "noise", "rep")
write.csv(config, "config.txt", row.names = FALSE)
