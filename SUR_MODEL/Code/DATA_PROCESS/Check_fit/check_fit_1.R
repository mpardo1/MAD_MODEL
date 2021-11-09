#DATA PROCESSING####
rm(list = ls())
library(gdata) 
library(ggplot2)
library(numbers)
library(tidyverse)
library(data.table)
library(multiplex)
library(tidyverse)
library("readxl")
library(reshape)
library(viridis)
library("ggpubr")

Path = "~/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_DETERMINISTIC/Output_Integration_final_fit.data"
int_sol_2 = data.frame(t(read.table(Path, header=FALSE)))

#### UPLOAD PARTICIPATION FILE #####
# Data of participation.
# Spain DATA:
Path_ages = "~/MAD_MODEL/SUR_MODEL/data/ages_days.csv"
ages = read.csv(Path_ages)
# Convert to data type date the registration time.
ages$date = as.Date(ages$date,"%Y-%m-%d") 
# Remove the data after the update of the app.
ages <- ages %>% filter(date < as.Date("2020-10-01","%Y-%m-%d") )
ages$time = as.numeric(ages$date - as.Date("2014-06-14","%Y-%m-%d") , units="days") + 1
ages$age_days = ages$age_days + 1

####UPLOAD REGISTRATION FILE ######
# Data with registration date with hour and minute and registration ID.
Path_reg = "~/MAD_MODEL/SUR_MODEL/data/ages_days_bcn.csv"
registration = read.csv(Path_reg)
registration <- registration %>% filter( registration$age_days == 0)
registration$date_reg <- as.Date(registration$date,'%Y-%m-%d') 
registration$age_days <- NULL
ref_date = min(registration$date)
# Convert to data type date the registration time.
# Remove the data after the update of the app.
init_date = "2014-06-14"
end_date = max(registration$date_reg)
registration <- registration %>% filter(date_reg < as.Date(end_date,"%Y-%m-%d") )
registration <- registration %>% filter(date_reg >= as.Date(init_date,"%Y-%m-%d") )
registration$boolean = 1
# Do a group by by the day of registration.
reg_group <- registration %>%  group_by(date_reg) %>% summarise(mean = mean(boolean), sum = sum(boolean), n = n())
reg_group_sort <- reg_group[order(reg_group$date_reg),]
# Erase the dummy columns.
reg_group_sort$mean <- NULL
reg_group_sort$sum <- NULL
reg_group_sort$diff <- NULL
reg_group_sort$time = as.numeric(reg_group_sort$date_reg - as.Date("2014-06-14","%Y-%m-%d") , units="days")
reg_group_sort <- reg_group_sort[,c(3,2)]
df_test = data.frame(time = c(0:max(reg_group_sort$time)), n =0)
# Añado los registros de los dias que no tengo para que tenga estos datos sin cortes.
reg_group_sort <- merge(x = df_test, y = reg_group_sort, by = "time", all.x = TRUE)
reg_group_sort$n.y[is.na(reg_group_sort$n.y)] <- 0
reg_group_sort$n = reg_group_sort$n.y
reg_group_sort$n.x <- NULL
reg_group_sort$n.y <- NULL
# Save FILE with DOWNLOADS DATA, A(t) en el modelo.
# write.dat(reg_group_sort, "/home/marta/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_DETERMINISTIC/")


###### PROCESS PARTICIPATION.DAT ########
# Create a matrix with the participant data in each row each age group dynamics.
ages_ord <- ages[with(ages, order(date, age_days)), ]
# Insert missing age_groups per date. 
# First create a column with time from 1 to end.
vec_uni = unique(ages_ord$date)
l_uni = length(vec_uni)-1
date_uni = data.frame(date = vec_uni, time =c(0:l_uni))
# Join the two data frames.
# merge_df <- merge(ages_ord,date_uni, by="date")
# dim_df = length(ages$date)
# merge_dt <- data.table(merge_df)
dim = max(ages_ord$time)
mat_ages <- matrix(0, dim+1, dim+1)
dim_l = length(ages_ord$date)
for(i in c(1:dim_l)){
  print(paste0("Index",i))
  mat_ages[ages_ord$age_days[i]+1,ages_ord$time[i]+1]=ages_ord$N[i]
}
mat_ages[1,] = c(0:dim)
length = length(mat_ages[1,])-2
# ###### SAVE FILE OBESERVED DATA-#########
# write.dat(mat_ages[c(1:length),c(1:length)], "/home/marta/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_CBL_ESTIMATION")
# write.dat(mat_ages, "/home/marta/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_CBL_ESTIMATION")

######  STANDARD DEVIATION #######

# All ages:
l = length(mat_ages[1,])
vec_sd <- seq(1,l,1)
vec_ages <- seq(1,l,1)
for(i in c(1:l)){
  vec_sd[i] = log10(sd(mat_ages[i,]))
}

df_sd <- data.frame(ages = vec_ages, sd = vec_sd)
ggplot(df_sd) + 
  geom_line(aes(ages,sd)) +
  xlab("Age of the participants set") + 
  ylab("Standard deviation") +
  theme_bw()

# All ages from the second one (fist one removed).
l = length(mat_ages[1,])
vec_sd <- seq(2,l,1)
vec_ages <- seq(2,l,1)
for(i in c(1:l-1)){
  vec_sd[i] = sd(mat_ages[i+1,])
}

df_sd <- data.frame(ages = vec_ages, sd = vec_sd)
ggplot(df_sd) + 
  geom_line(aes(ages,sd)) +
  xlab("Age of the participants set") + 
  ylab("Standard deviation")+
  theme_bw()


mat_diff <- mat_ages[2,] - mat_ages[1,]
df_sd <- data.frame(ages = vec_ages, sd = vec_sd)
ggplot(df_sd) + 
  geom_line(aes(ages,sd)) +
  xlab("Age of the participants set") + 
  ylab("Standard deviation")+
  theme_bw()



#----------------------------------------------------------------------#
mat_t = data.frame(t(mat_ages[c(1,2),c(1:30)]))
len = length(mat_ages[1,])-1
mat_ages_1 = rbind(c(0:len),mat_ages)
vec_i = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,
          24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40,
          41,42,43,44,45,46,47,48,49,50,51,52,53,54,
          55,56,57,58,59,60,1730,1731,1732,1733,1734,1735,1736,1737,1738,1739,1740,1741,1742,1743,1744,
          1745,1746,1747,1748,1749,1780,1781,1782,1783,1784,1785,1786,1787,1788,1789,1790)
mat_fil = data.frame(t(mat_ages)[,vec_i])

rm(ages)
rm(ages_ord)
rm(date_uni)
rm(mat_ages_1)

##### PLOTS TEMPORAL SERIES PARTICIPATION ########
###### UPLOAD INTEGRATION FILES ######
# Data of participation.
Path = "~/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_DETERMINISTIC/Output_Integration_final_fit.data"
# Path = "/home/marta/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_DETERMINISTIC/Output_Integration.dat"
# 
# Mac PATH
#Path = "/Users/celsaaraujobarja/Documents/PhD/Output_Int/Output_Integration.dat"
int_sol = data.frame(t(read.table(Path, header=FALSE)))

#################COMPUTE ERROR#############################
max_age <- length(int_sol[1,])
error <- matrix(2,max_age,2)
max_time <- min(max(mat_ages[1,]), max(int_sol[,1]))

mat_err <- matrix(1, max_age,max_time)
for(i in c(1:max_age-1)){
  print(paste0("i:",i))
  mat_err[i,] <- abs(mat_ages[i+1,1:max_time] - int_sol[1:max_time,i+1])
  error[i,1] <- mean( mat_err[i,])
  error[i,2] <- sd( mat_err[i,])
}

hist(mat_err[56,], breaks=12, col="red") 

mean_ages_err <- seq(1,max_time,1)
vec <- c(1:max_time)
for(i in c(1:max_time)){
  mean_ages_err[i] <- mean(mat_err[,i])
}

df_err_mean <- data.frame(time = vec, mean_error_tot <- mean_ages_err)
ggplot(df_err_mean) + 
  geom_line(aes(vec,mean_ages_err))

df_err <- data.frame(ages = c(1:max_age), error_mean= error[,1])
ggplot(df_err) + 
  geom_line(aes(ages,error_mean))

df_err <- data.frame(ages = c(1:max_age), error_mean= error[,2])
ggplot(df_err) + 
  geom_line(aes(ages,error_mean))

hist(df_err$error_mean, breaks=12, col="red") 

# int_sol = data.frame(int_sol[,vec_i])
int_sum = rowSums(int_sol[1:2300,2:length(int_sol[1,])])
mat_sum = rowSums(mat_fil[1:2300,2:length(mat_fil[1,])-1])
df_sum <- data.frame(Time = int_sol[1:2300,1], Simulation = int_sum ,Observed = mat_sum)
df_sum <- reshape2::melt(df_sum, id.vars = c("Time"))

# Plot number of participants at each time.
ggplot(df_sum,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  ylab("Number of participants") +
  scale_color_manual(values=c('#9E329F','#1642FE'))+
  theme(text = element_text(size=20))

# Reshape dataframe for ggplot:
df_join <- merge(int_sol, mat_fil, by = "X1")

# First Age group:
theme_set(
  theme_bw() +
    theme(text = element_text(size=12)) + 
    theme(legend.title = element_blank())
)

df_aux <- df_join[,c("X1","X2.y","X2.x")]
colnames(df_aux) <- c("Time","Observed","Simulation")
df_plot <- reshape2::melt(df_aux, id.vars = c("Time"))

# 1 day old:
PLOT_1 <- ggplot(df_plot,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F','#1642FE'))+
  ylab("Number of participants")+
  xlab("Time (days)") 
PLOT_1
df_plot <- df_plot %>% filter(variable == "Observed" )
ggplot(df_plot,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F'))+
  ylab("Number of participants")+
  xlab("Time (days)") 

rm(df_plot)
rm(df_aux)
df_aux <- df_join[,c("X1","X26.y","X26.x")]
colnames(df_aux) <- c("Time","Observed","Simulation")
df_plot <- reshape2::melt(df_aux, id.vars = c("Time"))

# 25 day old:
PLOT_2 <- ggplot(df_plot,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F','#1642FE'))+
  ylab("Number of participants")+
  xlab("Time (days)") 
PLOT_2
df_plot <- df_plot %>% filter(variable == "Observed" )

ggplot(df_plot,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F'))+
  ylab("Number of participants")+
  xlab("Time (days)")

# Second Age group:
df_aux <- df_join[,c("X1","X36.y","X36.x")]
colnames(df_aux) <- c("Time","Observed","Simulation")
df_plot <- reshape2::melt(df_aux, id.vars = c("Time"))

# 36 day old:
PLOT_3 <- ggplot(df_plot,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F','#1642FE'))+
  ylab("Number of participants")+
  xlab("Time (days)") 
PLOT_3
df_plot <- df_plot %>% filter(variable == "Observed" )
ggplot(df_plot,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F'))+
  ylab("Number of participants")+
  xlab("Time (days)") 

rm(df_plot)
rm(df_aux)
df_aux <- df_join[,c("X1","X61.y","X61.x")]
colnames(df_aux) <- c("Time","Observed","Simulation")
df_plot <- reshape2::melt(df_aux, id.vars = c("Time"))

# 61 day old:
ggplot(df_plot,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F','#1642FE'))+
  ylab("Number of participants")+
  xlab("Time (days)") 

# Third Age group:
df_aux <- df_join[,c("X1","X62.y","X62.x")]
colnames(df_aux) <- c("Time","Observed","Simulation")
df_plot <- reshape2::melt(df_aux, id.vars = c("Time"))

# 1730 day old:
PLOT_4 <- ggplot(df_plot,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F','#1642FE'))+
  ylab("Number of participants")+
  xlab("Time (days)")
PLOT_4
df_plot <- df_plot %>% filter(variable == "Observed" )
ggplot(df_plot,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F','#1642FE'))+
  ylab("Number of participants")+
  xlab("Time (days)")

rm(df_plot)
rm(df_aux)
df_aux <- df_join[,c("X1","X81.y","X81.x")]
colnames(df_aux) <- c("Time","Observed","Simulation")
df_plot <- reshape2::melt(df_aux, id.vars = c("Time"))

# 1748 day old:
ggplot(df_plot,aes(Time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F','#1642FE'))+
  ylab("Number of participants")+
  xlab("Time (days)")

figure <- ggarrange(PLOT_1, PLOT_2, PLOT_3, PLOT_4,
                    common.legend = TRUE, legend = "top",
                    ncol = 2, nrow = 2)
figure
###### PLOTS AGE DISTRIBUTION ######

#mat_fil_t = data.frame(t(mat_fil[c(1:31),c(2:32)]))
#mat_fil_t$time = c(0:30)
init_sol_mat = array(int_sol)
# DF to check wheter the downloads are the participants of age 0.
check_df = data.frame(time = mat_ages[1,1:31], P_0 = mat_ages[2,1:31], down = reg_group_sort$n[1:31] )
check_df$diff = check_df$down - check_df$P_0
int_sol_t = data.frame(t(int_sol[c(1:31),c(6:36)]))
int_sol_t$time = c(0:30)
##Check whether the age distribution is correct.
time = 30
mat_fil_t = data.frame(time = mat_ages[1,1:31], Ob_30 = mat_ages[2:32,time], te_30 = t(int_sol[time,6:36]))
ggplot(mat_fil_t) +
  geom_line(aes(time,Ob_30), color = "red") +
  geom_line(aes(time,V30),size = 0.6) +
 # xlim(1756,2000) + 
  xlab("Age of the participant") + 
  ylab("Number of participants") 
