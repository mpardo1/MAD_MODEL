rm(list = ls())
library(gdata) 
library(segmented)
library(e1071)
library(ggplot2)


################ READ DATA ######################

Path = "/home/marta/Documentos/PHD/2021/SUR_Model/PARAM_CONFIG/"
Path = "/home/marta/Documentos/PHD/2021/SUR_Model/RESULTS_ESTIMATION/UNTIL_TODAY/"
setwd(Path)
file_list <- list.files(Path)
file_list <- file_list[names(file_list) != "NOT_RAND" && names(file_list) != "RAND" && names(file_list) != "UNTIL_TODAY"]
for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("dataset")){
    print("First choice")
    print(file)
    dataset <- read.table(file, header=TRUE, sep="\t")
    file_list <- file_list[names(file_list) != file]
  }else{
    print("Second choice")
    print(file)
    temp_dataset <- read.table(file, header=TRUE)
    dataset <- rbind(dataset, temp_dataset)
    rm(temp_dataset)
  }
  file_list <- file_list[names(file_list) != file]  
}
################ PROCESS DATA ######################
#dataset <- read.table(Path, header = TRUE)
dataset <- dataset[order(dataset$NegLogLike),]
write.csv(dataset, "ORD_CONF_PARAM")

Path = "~/Documentos/RESULTS_ESTIMATION/NOT_RAND"
Path = "/home/marta/Documentos/Output_Int/PDE"
