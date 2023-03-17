# setwd("C:/Users/athen/Research/Codes/repo_PHD/Inference_in_natural_sciences/Inverse_Problems/")
setwd("/storage/homefs/ag19v929/Research/Codes/repo_PHD/Inference_in_natural_sciences/Inverse_Problems")

list_files <- list.files(path = "./res/scores/")
res <- data.frame()
res_vanilla <- data.frame()

for(ifile in seq_along(list_files)){
  file <- list_files[ifile]
  
  if(ifile%%50==1){
    cat(round(100*ifile/length(list_files), 1), "%\n")
  }
  temp <- read.csv(paste0("./res/scores/", file))
  temp$model <- strsplit(file, "_sampling_")[[1]][1]
  temp$strategy <- strsplit(strsplit(file, "_sampling_")[[1]][2], "_jobnumber_")[[1]][1]
  temp$X <- NULL
  if(ncol(temp)==408){
    res_vanilla <- rbind(res_vanilla, temp)
  }else{
    res <- rbind(res, temp)
  }
}
write.csv(res, file="./res/combined_scores.txt", row.names = FALSE)
write.csv(res_vanilla, file="./res/combined_scores_vanilla.txt", row.names = FALSE)