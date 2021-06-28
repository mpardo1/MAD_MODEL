rm(list = ls())
library(easypackages)
libraries("gdata", "ggplot2", "numbers","tidyverse","data.table","multiplex","reshape","viridis","stats","ggpubr","ggstatsplot","e1071","mlr3misc","deSolve", "gganimate") 


# Different path for MAc and Ubuntu.
PC = "/Users/celsaaraujobarja/Documents"
#PC = "/home/marta/Documentos"

####INPUT#####
#Temperatures from Barcelona :
Path_temp = "/PHD/2021/Mosquito_model/data/bcn_weather_daily.Rds"
Path_temp = paste(PC,Path_temp, sep="")

temp <-read_rds(Path_temp)
temp$date = as.Date(temp$date , "%Y-%m-%d")
temp <- temp %>%  group_by(date) %>% summarise(mean_temp = mean(valor))

ggplot(temp) +
  geom_line(aes(date, mean_temp))+
  ggtitle("Mean temperature Barcelona") + 
  xlab("Mean temperature")

#######RHO(M)#######
#Path= "/home/marta/Documentos/PHD/2021/SUR_Model/Code/OUTPUT/df_rho.dat"
Path_rel= "/PHD/2021/SUR_Model/Code/OUTPUT/df_rho.dat"
Path = paste(PC,Path_rel, sep="")
df_rho <- data.frame(t(read.table(Path, header=FALSE)))
colnames(df_rho) <- c("time", "rho", "date")
df_rho$date = as.Date(df_rho$date , "%Y-%m-%d")
df_rho$time <- NULL
df_rho$rho <- as.numeric(df_rho$rho)
ggplot(df_rho) + 
  geom_line(aes(x = date, y =rho)) +
  ggtitle("Encounter rate computed from Citizen Science model") 

####  Functional responses####
# Gonotrophic cycle
gonot <- function(T){
  a = 0.045*T^2-2.617*T+41.105
  val = min(1,1/a)
  return(val)
}

# Development rate
d_L <- function(T){
  val = 1/(0.14457*T^2 - 8.24857*T + 124.80857)
  return(val)
}

# Larva mortality rate
delta_L <- function(T){
  a = abs(-0.1305*T^2+3.868*T+30.83)
  val = min(1, 1/a)
  return(val)
}

# Adult mosquito mortality rate
delta_A <- function(T){
  a = abs(-0.1921*T^2+8.147*T-22.98)
  val =  min(1, 1/a)
  return(val)
}



vec = seq(0,40,1)
df_gonot_vec <- data.frame(temp = vec, gonot =unlist(lapply(vec,gonot)))
df_dL_vec <- data.frame(temp = vec, dL = unlist(lapply(vec,d_L)))
df_deltaL_vec <- data.frame(temp = vec, deltaL = unlist(lapply(vec,delta_L)))
df_deltaA_vec <- data.frame(temp = vec, deltaA = unlist(lapply(vec,delta_A)))

