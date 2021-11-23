# DATA PROCESSING####
rm(list = ls())
library(ggplot2)
library(numbers)
library(tidyverse)
library(deSolve)

# UPLOAD FILES ####
# Data of participation.
Path_ages = "~/MAD_MODEL/SUR_MODEL/data/ages_days.csv"
ages = read.csv(Path_ages)
# Convert to data type date the registration time.
ages$date = as.Date(ages$date,"%Y-%m-%d") 
length_reg = length(ages$date)

# REGSITRATION DATA ###
registration <- read.csv(Path_ages)
registration <- registration %>% filter( registration$age_days == 0)
# Convert to data type date the registration time.
registration$date <- as.Date(registration$date,'%Y-%m-%d') 
min_date <- min(registration$date)
max_date <- max(registration$date)
df_data <- data.frame(date = seq(min_date, max_date, "days"))
reg_df <- merge(df_data,registration, by="date", all = TRUE)
reg_df[is.na(reg_df)] <- 0
reg_df$age_days <- NULL
for(i in c(1:(length(registration$date)-1))){
  a <- as.numeric(registration$date[i+1] - registration$date[i])
  if(a > 1){
    print(paste0("diff date", registration$date[i]))
  }
}


reg_df$time <- as.numeric(reg_df$date - min_date)
registration_df <- t(as.matrix(data.frame( time = reg_df$time, reg_df$N)))
Path <- paste0( "~/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_CBL_ESTIMATION/Downloads_", Sys.Date() ,"_", max(reg_df$time),".dat")
write.table(registration_df, Path,sep="\t",col.names = FALSE,row.names = FALSE)
#-----------------------------------------------------------------------------------#
## Create the matrix with observed data, each row one group class each column each time.
min_date <- min(ages$date)
ages$time <- as.numeric(ages$date - min_date)
dim <- max(ages$time) + 1
mat_ages <- matrix(0, dim, dim)
d <- nrow(ages)
for(i in c(1:d)){
  mat_ages[ages$age_days[i],ages$time[i]]=ages$N[i]
}
mat_ages <- rbind(seq(1,ncol(mat_ages),1),mat_ages)
mat_ages <- cbind(replicate(nrow(mat_ages), 0), mat_ages)
# Save matrix to do the estimation parameter problem in C:
Path <- paste0( "~/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_CBL_ESTIMATION/Observed_", Sys.Date() ,"_", ncol(mat_ages),".dat")
write.table(mat_ages, Path,sep="\t",col.names = FALSE,row.names = FALSE)

Path <- paste0( "~/MAD_MODEL/SUR_MODEL/data/Observed_", Sys.Date() ,"_", ncol(mat_ages),".dat")
write.table(mat_ages, Path,sep="\t",col.names = FALSE,row.names = FALSE)
# Create a data frame with the observed data with the time and date:
mat_ages <- t(mat_ages)
min_date <- min(ages$date) - 1
max_date <- max(ages$date)
ob_data <- as.data.frame(mat_ages)
date <- seq(min_date,max_date,"days")
ob_data <- cbind(date, ob_data)

Path = "~/MAD_MODEL/SUR_MODEL/data/Ob_data.rds"
saveRDS(ob_data,Path)

#------------------------------------------------------------------------------------#
#### Observed data for Barcelona #####
Path_ages = "~/MAD_MODEL/SUR_MODEL/data/ages_days_bcn.csv"
ages = read.csv(Path_ages)
# Convert to data type date the registration time.
ages$date = as.Date(ages$date,"%Y-%m-%d") 
length_reg = length(ages$date)
## Create the matrix with observed data, each row one group class each column each time.
min_date <- min(ages$date)
ages$time <- as.numeric(ages$date - min_date)
dim <- max(ages$time) + 1
mat_ages <- matrix(0, dim, dim)
d <- nrow(ages)
for(i in c(1:d)){
  mat_ages[ages$age_days[i],ages$time[i]]=ages$N[i]
}
mat_ages <- rbind(seq(1,ncol(mat_ages),1),mat_ages)
mat_ages <- cbind(replicate(nrow(mat_ages), 0), mat_ages)

mat_ages <- t(mat_ages)
min_date <- min(ages$date) - 1
max_date <- max(ages$date)
ob_data <- as.data.frame(mat_ages)
date <- seq(min_date,max_date,"days")
ob_data <- cbind(date, ob_data)

Path = "~/MAD_MODEL/SUR_MODEL/data/Ob_data_bcn.rds"
saveRDS(ob_data,Path)

#-----------------------------------------------------------------------------------#
sum <- array(0,dim)
for(i in c(1:dim)){
  sum[i] = sum(mat_ages[i,])
  if(sum[i] == 0){
    print(paste0("Index",paste0( i, "with 0 all")))
  }
}
write.dat(mat_ages[1:4,1:12], "/home/marta/Documentos/SUR Model/Code/RStudio/Fitting/data/mat_ages_12P_4R.dat")
len = length(mat_ages[1,])-1
mat_ages_1 = rbind(c(0:len),mat_ages)
write.dat(mat_ages_1[1:60,1:12], "/home/marta/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_CBL_ESTIMATION/mat_ages_12P_4R_1.dat")
vec_i = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30,31,81,131,181,231,281,331,381,431,481,531,581,631,681,731,781,831,881,931,981,1206,1316,1426,1536,1646,1756,1866,1976,2086,2196)
write.dat(mat_ages_1[1:2378,1:2378], "/home/marta/Documentos/mat_ages_12P_4R_1.dat")
mat_fil = mat_ages_1[vec_i,1:2378]
remove(ages)
remove(date_uni)
# remove(mat_ages)
# remove(mat_ages_1)
remove(merge_df)
remove(reg_group_sort)
remove(merge_dt)

# -----------------PLOTS--------------------#####
# Create a dataframe with the filter. To do the animated plot.
df_filt <- ages_ord[ which( ages_ord$date > "2016-01-01" & ages_ord$date < "2017-01-01" & ages_ord$age_days < 4), ]
plot <- ggplot(df_filt) +
  geom_line(aes(x = age_days, y = N),size = 0.5) +
  theme(text = element_text(size=20)) +
  ylim(0,500) +
  # scale_x_continuous(breaks = c(2,3,4),
  #                    labels = c("$400","$4 000","$40 000"),
  #                    limits = c(2,4.1)) +
  transition_states(date, transition_length = 1, state_length = 1) +
  ylab('Number of participants') +
  xlab('Participant age') +
  ggtitle('Day: {closest_state}')
animate(plot, render = gifski_renderer() ,width =400, height = 400, nframes = 480, fps = 24)

anim1 <- ggplot(df_filt, aes(age_days, N)) +
  geom_point() +
  labs(title = "{closest_state}") +
  transition_states(date, transition_length = 1, state_length = 1) +
  scale_y_continuous(limits = c(0,100))+
  enter_fade() +
  exit_fade()

max_f = 100
df_file <- data.frame(matrix(0,ncol = 3, nrow = max_f),col.type = c("integer","Date","numeric"))
df_file$X2 = "2020-01-01"
df_file %>% mutate(X2 = as.Date(X2,"%Y-%m-%d"))
for(i in c(1:max_f)){
  df_filt1 <- ages_ord[ which( ages$date == date_uni$date[i]), ]
  ggsave(paste0("/home/marta/Documentos/SUR Model/Code/RStudio/Fitting/PLOT_DIST_AGES/plot",paste0(df_filt1$date[i],".png")))
    df_file$X1[i] = i
    df_file$X2[i] = as.Date(df_filt1$date[i],"%Y-%m-%d")
    df_file$X3[i] = paste0("plot",paste0(df_filt1$date[i],".png"))
  plot1 <- ggplot(df_filt1) + 
    geom_line(aes(x = age_days, y = N),size = 0.5) 
}

locations <- unique(df_file$X1)

for(i in 1:length(locations)) {
  images <- map(df$file[df$loc == locations[i]], image_read)
  images <- image_join(images)
  animation <- image_animate(images, fps = 1)
  image_write(animation, paste0(locations[i], ".gif"))
}

list.files(path='/home/marta/Documentos/SUR Model/Code/RStudio/Fitting/PLOT_DIST_AGES/', pattern = '*.png', full.names = TRUE) %>% 
  image_read() %>% # reads each path file
  image_join() %>% # joins image
  image_animate(fps=4) %>% # animates, can opt for number of loops
  magick::image_write("FileName.gif") # write to current dir

df_filt1 <- ages_ord[ which( ages$date == "2020-01-01"), ]
plot1 <- ggplot(df_filt1) + 
  geom_line(aes(x = age_days, y = N),size = 0.5) 
# 
# # # Another way of doing an animated plot.
# p <- df_filt %>%
#   plot_ly(
#     x = ~age_days,
#     y = ~N,
#     frame = ~date,
#     hoverinfo = "text",
#     type = 'scatter',
#     mode = 'markers'
#   ) %>%
#   layout(
#     xaxis = list(
#       type = "log"
#     )
#   )
# png("myPlot.png")
# plot(rnorm(1000),rnorm(1000))
# dev.off()
# library(png)
# P1 <- readPNG("myPlot.png")
# library(caTools)
# write.gif(P1,"myPlot.gif")
# showGIF <- function(fn) system(paste("display",fn))
# showGIF("myPlot.gif")
# unlink("myPlot.gif") 
