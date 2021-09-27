rm(list = ls())
library(easypackages)
libraries("gdata", "ggExtra","ggplot2", "numbers","tidyverse",
          "data.table","multiplex","reshape","viridis","stats",
          "ggpubr","ggstatsplot","e1071","mlr3misc","deSolve", "gganimate") 

# Moving average function:
mov_avg <- function(dim, dt, vec_in){
  vec_out = c(1:dim)
  for(i in c(1:dim)){
    if( i < (dim - dt)){
      end = i + dt
      vec_out[i] = mean(vec_in[i:end])
    }else{
      vec_out[i] = mean(vec_in[i:dim])
    }
  }
  return(vec_out)
}

###--------------------BG COUNTS------------------------####
Path = '~/MAD_MODEL/MAD_MODEL/data/bcn_bgcount_time_profile.csv'

bg_traps <- read.csv(Path)
bg_traps$date <- as.Date(bg_traps$date,'%Y-%m-%d') 
ggplot(bg_traps, aes(x = date, y = value)) +
  geom_point()
pop_density = 1664182
bg_traps$prop <- bg_traps$value/pop_density
ggplot(bg_traps) + 
  geom_line(aes(date, density))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab("Number of mosquito adults per human (BG estimation)") +
  theme_bw()

l = length(bg_traps$value)
dt = 7
mov = mov_avg(l,dt,bg_traps$value)
df_bg_mov <- data.frame(date = bg_traps$date, bg_counts = mov)
ggplot(df_bg_mov) + 
  geom_line(aes(date, bg_counts))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab("Number of mosquito adults per human (BG estimation)") +
  theme_bw()

#####------------------------SCATTER PLOT---------------------------#######
Path = '~/MAD_MODEL/SUR_MODEL/Code/adults.rds'
adults <- readRDS(Path)
ggplot(adults) + 
  geom_line(aes(date, A))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab("Number of mosquito adults per human (BG estimation)") +
  theme_bw()
l = length(adults$A)
adults$mov_A <- mov_avg(l,dt,adults$A)
df_tot <- merge(adults, df_bg_mov, by = "date")

# Scatter plot with marginal distributions for all data.
ggscatterstats(data = df_tot,
               x = bg_counts,
               y = mov_A, 
               xlab ="Adult mosquito  (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

# Scatter plot by seasons:
df_season1 <- df_tot %>% filter(date > as.Date("2018-05-15","%Y-%m-%d") & date < as.Date("2018-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season1,
               x = bg_counts,
               y = mov_A, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#Season 2:
df_season2 <- df_tot %>% filter(date > as.Date("2019-04-01","%Y-%m-%d") & date < as.Date("2019-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season2,
               x = bg_counts,
               y = mov_A, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#Season 3:
df_season3 <- df_tot %>% filter(date > as.Date("2020-04-01","%Y-%m-%d") & date < as.Date("2020-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season3,
               x = bg_counts,
               y = mov_A, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

