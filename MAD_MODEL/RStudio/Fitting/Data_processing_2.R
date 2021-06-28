# DATA PROCESSING####
rm(list = ls())
library(gdata) 
library(segmented)
library(e1071)
library(ggplot2)
library(optimx)
library(zoom)
library(mlr3misc)
library(numbers)
library(tidyverse)
library(SparkR)
library(deSolve)
library(coda)
library(rootSolve)
library(FME)
library(data.table)
library(multiplex)
library(gganimate)
library(ggflags)
library(gifski)
library(multiplex)
library(png)
library(caTools)
library(plotly)
library(gapminder)
library(purrr)
library(magick)
library(magrittr)

# UPLOAD FILES ####
# Data of participation.
Path_ages = "/home/marta/Documentos/SUR Model/Code/RStudio/Fitting/data/ages_days.csv"
ages = read.csv(Path_ages)
# Convert to data type date the registration time.
ages$date = as.Date(ages$date,"%Y-%m-%d") 
length_reg = length(ages$date)

# Data with registration date with hour and minute and registration ID.
Path_users = "/home/marta/Documentos/SUR Model/Code/RStudio/Fitting/data/register_data_tigausers.csv"
registration = read.csv(Path_users)

# Downloads .data : Downloads_Transposed_2378
Path_down = "/home/marta/Documentos/SUR Model/Code/RStudio/Fitting/data/Downloads_Transposed_2378.data"
down_process = t(read.table(Path_down))
l_down = max(down_process[,1])
down_process_full = matrix(0, l_down, 2)
down_process_full[,1]
for(i in c(1:l_down)){
  down_process_full[i,1]
}
# Convert to data type date the registration time.
registration$date_reg = as.Date(registration$registration_time,'%Y-%m-%d %H:%M:%S') 
registration$boolean = 1
# Do a group by by the day of registration.
reg_group <- registration %>%  group_by(date_reg) %>% summarise(mean = mean(boolean), sum = sum(boolean), n = n())
reg_group_sort <- reg_group[order(reg_group$date_reg),]
# Erase the dummy columns.
reg_group_sort$mean <- NULL
reg_group_sort$sum <- NULL
reg_group_sort$diff <- NULL
len = length(reg_group_sort$n)-1
reg_group_sort$time = 0
for(i in c(1:len)){
  diff  = reg_group_sort$date_reg[i+1] - reg_group_sort$date_reg[i] 
  reg_group_sort$time[i+1] = reg_group_sort$time[i] + diff
  print(paste0("Date difference:",diff))
  print(paste0("Time i+1:",reg_group_sort$time[i+1]))
}
reg_group_sort$date_reg <- NULL
reg_group_sort <- t(reg_group_sort)
reg_group_sort_1 <- matrix(0,2,len+1)
reg_group_sort_1[1,] = reg_group_sort[2,]
reg_group_sort_1[2,] = reg_group_sort[1,]
write.table(reg_group_sort_1, "~/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_CBL_ESTIMATION/Downloads_2353_test.dat",sep="\t",col.names = FALSE,row.names = FALSE)
# Remove all data frames.
remove(reg_group)
remove(registration)

##############-----------------------PROCESS .DAT-----------------------###########
# Create a matrix with the participant data in each row each age group dynamics.
ages_ord <- ages[with(ages, order(date, age_days)), ]
# Insert missing age_groups per date. 
# First create a column with time from 1 to end.
vec_uni = unique(ages_ord$date)
l_uni = length(vec_uni)-1
date_uni = data.frame(date = vec_uni, time =c(0:l_uni))
# Join the two data frames.
merge_df <- merge(ages_ord,date_uni, by="date")
dim_df = length(ages$date)
merge_dt <- data.table(merge_df)
dim = max(merge_dt$time)
mat_ages <- matrix(0, dim, dim)
dim_l = length(merge_dt$date)
for(i in c(1:dim_l)){
  print(paste0("Index",i))
  mat_ages[merge_dt$age_days[i],merge_dt$time[i]]=merge_dt$N[i]
}
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
