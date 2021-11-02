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
library(png)
library(caTools)

# Data of participation.
Path_ages = "/home/marta/Documentos/SUR Model/Code/RStudio/Fitting/data/ages_days.csv"
ages = read.csv(Path_ages)
# Convert to data type date the registration time.
ages$date = as.Date(ages$date,"%Y-%m-%d") 
length_reg = length(ages$date)

# Data with registration date with hour and minute and registration ID.
Path_users = "/home/marta/Documentos/SUR Model/Code/RStudio/Fitting/data/register_data_tigausers.csv"
registration = read.csv(Path_users)
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
# Remove all data frames.
remove(reg_group)
remove(registration)

##############-----------------------PROCESS .DAT-----------------------###########
# Create a matrix with the participant data in each row each age group dynamics.
ages_ord <- ages[with(ages, order(date, age_days)), ]
# Insert missing age_groups per date. 
# First create a column with time from 1 to end.
vec_uni = unique(ages_ord$date)
l_uni = length(vec_uni)
date_uni = data.frame(date = vec_uni, time =c(1:l_uni))
# Join the two data frames.
merge_df <- merge(ages_ord,date_uni, by="date")
dim_df = length(ages$date)
min_i = 1
max_i = dim_df-1
for(i in c(min_i:max_i)){
  diff_age = as.numeric(merge_df$age_days[i+1]) - as.numeric(merge_df$age_days[i])
  if(merge_df$date[i] == merge_df$date[i+1] && (diff_age > 1)){
    print(paste0("Gap between the ages at i:", i))
    date_i = merge_df$date[i]
    days_i1 = as.numeric(merge_df$age_days[i])
    time_i = merge_df$time[i]
    for(j in c(1:diff_age)){
      if(is.na(date_i) | is.na(days_i1)){
        print("******************NA FOUND*******************")
        exit()
      }
      merge_df <- rbind(merge_df,c(as.character(date_i),days_i1+j,0,time_i))
    }
  }
}
merge_df <- ages[with(merge_df, order(date, age_days)), ]
merge_dt <- as.data.frame.table(merge_df)
dim = max(ages$age_days)
vec = ages$N
mat_ages <- matrix(0, dim, dim)
ind = 1
for(i in c(0:dim)){
  for(j in c(0:dim)){
    if(j<=i){
      # print("index i:")
      # print(i)
      print(" \n")
      print(paste0("Index:",i,",",j))
      print(paste0("Index vec:",ind+j-1))
      print(paste0("Value vec:",vec[ind+j-1]))
      mat_ages[i,j] = vec[ind+j-1]
      if(j==i){ind = ind+j}
    }else{
      # print("0")
      mat_ages[i,j] = 0
    }
  }
}

# date <- unique(ages_ord$date_1)
# l_vec = length(vec_dates)
# time <- c(1:l_vec)
# df_ages <- data.frame(date,time)
##############--------------------------PLOTS-----------------------###########
# Create a new data frame, applying some filters in the ages df.
df_filt <- ages[ which( (ages$age_days == 2 | ages$age_days == 150 | ages$age_days ==  250) & ages$date > "2016-01-01" & ages$date < "2017-01-01") , ]
down_filt <- reg_group_sort[ which(reg_group_sort$date_reg> "2016-01-01"  & reg_group_sort$date_reg < "2017-01-01"),]
# Change the column age_days to string to be understood as a category.
df_filt$age_days = as.character(df_filt$age_days)

# Make a plot with the downloads and participation dynamics
# plot_ages <- ggplot() +
#   geom_line(aes(x = date, y= N, color=age_days),df_filt, size=0.4)+
#   geom_line( aes(x=date_reg, y=n, color = "Downloads"), down_filt, size = 0.4) +
#  # scale_y_continuous(name="Number of participants") + 
#   ggtitle("Distribution of participants along time")+
#   #labs(fill = "Age group") +
#   theme_bw()+
#   labs(x = "Date",
#        y = "Number of participants",
#        color = "Legend") +
#   scale_color_manual(values = c("blue","green","red","black"))
# plot_ages

# Remove old df.
remove(df_filt)
remove(down_filt)
remove(plot_ages)
remove(reg_group_sort)
