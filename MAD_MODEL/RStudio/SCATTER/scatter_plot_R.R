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
  ylab("Number of mosquito df_rho per human (BG estimation)") +
  theme_bw()

l = length(bg_traps$value)
dt = 7
mov = mov_avg(l,dt,bg_traps$value)
df_bg_mov <- data.frame(date = bg_traps$date, bg_counts = mov)
ggplot(df_bg_mov) + 
  geom_line(aes(date, bg_counts))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab("Number of mosquito df_rho per human (BG estimation)") +
  theme_bw()

#####------------------------SCATTER PLOT---------------------------#######

####### ADULTS SIMULATIONS ######
Path = '~/MAD_MODEL/SUR_MODEL/Code/adults_sim.rds'
vector_sim <- readRDS(Path)
ggplot(vector_sim) + 
  geom_line(aes(date, Adults))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab("Number of mosquito adults per human (BG estimation)") +
  theme_bw()
l = length(vector_sim$Adults)
vector_sim$mov_A <- mov_avg(l,dt,vector_sim$Adults)
df_tot <- merge(vector_sim, df_bg_mov, by = "date")

# Scatter plot with marginal distributions for all data.
ggscatterstats(data = df_tot,
               x = bg_counts,
               y = mov_A, 
               xlab ="Adult mosquito  (BG traps)" , 
               ylab = "Adult mosquito Simulations") 

ggscatterstats(data = df_tot,
               x = mov_A,
               y = bg_counts, 
               xlab ="Adult mosquito Simulations" , 
               ylab = "Adult mosquito  (BG traps) ") 

# Scatter plot by seasons:
df_season1 <- df_tot %>% filter(date > as.Date("2018-05-15","%Y-%m-%d") 
                                & date < as.Date("2018-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season1,
               x = bg_counts,
               y = mov_A, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#Season 2:
df_season2 <- df_tot %>% filter(date > as.Date("2019-04-01","%Y-%m-%d") 
                                & date < as.Date("2019-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season2,
               x = bg_counts,
               y = mov_A, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#Season 3:
df_season3 <- df_tot %>% filter(date > as.Date("2020-04-01","%Y-%m-%d") 
                                & date < as.Date("2020-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season3,
               x = bg_counts,
               y = mov_A, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#--------------------------------------------------------------------------------------------#
##### RHO SIMULATIONS ######
Path_rho <- "~/MAD_MODEL/SUR_MODEL/Code/rho_sim.rds"
df_rho <- readRDS(Path_rho)
ggplot(df_rho) + 
  geom_line(aes(date, rho))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab(expression(rho)) +
  theme_bw()
l = length(df_rho$rho)
df_rho$mov_rho <- mov_avg(l,dt,df_rho$rho)
df_tot <- merge(df_rho, df_bg_mov, by = "date")

# Scatter plot with marginal distributions for all data.
ggscatterstats(data = df_tot,
               x = bg_counts,
               y = mov_rho, 
               xlab ="Adult mosquito  (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

# Scatter plot by seasons:
df_season1 <- df_tot %>% filter(date > as.Date("2018-05-15","%Y-%m-%d") 
                                & date < as.Date("2018-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season1,
               x = bg_counts,
               y = mov_rho, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#Season 2:
df_season2 <- df_tot %>% filter(date > as.Date("2019-04-01","%Y-%m-%d") 
                                & date < as.Date("2019-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season2,
               x = bg_counts,
               y = mov_rho, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#Season 3:
df_season3 <- df_tot %>% filter(date > as.Date("2020-04-01","%Y-%m-%d") 
                                & date < as.Date("2020-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season3,
               x = bg_counts,
               y = mov_rho, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#--------------------------------------------------------------------------------------------#
##### RHO OBSERVED ######
Path_rho_ob <- "~/MAD_MODEL/SUR_MODEL/Code/rho_observed.rds"
df_rho_ob <- readRDS(Path_rho_ob)
ggplot(df_rho_ob) + 
  geom_line(aes(date, rho))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab(expression(rho)) +
  theme_bw()
l = length(df_rho_ob$rho)
df_rho_ob$mov_rho_1 <- mov_avg(l,dt,df_rho_ob$rho)
df_tot <- merge(df_rho_ob, df_bg_mov, by = "date")

# Scatter plot with marginal distributions for all data.
ggscatterstats(data = df_tot,
               x = bg_counts,
               y = mov_rho, 
               xlab ="Adult mosquito  (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

# Scatter plot by seasons:
df_season1 <- df_tot %>% filter(date > as.Date("2018-05-15","%Y-%m-%d") 
                                & date < as.Date("2018-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season1,
               x = bg_counts,
               y = mov_rho, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#Season 2:
df_season2 <- df_tot %>% filter(date > as.Date("2019-04-01","%Y-%m-%d") 
                                & date < as.Date("2019-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season2,
               x = bg_counts,
               y = mov_rho, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#Season 3:
df_season3 <- df_tot %>% filter(date > as.Date("2020-04-01","%Y-%m-%d") 
                                & date < as.Date("2020-11-01","%Y-%m-%d") )
ggscatterstats(data = df_season3,
               x = bg_counts,
               y = mov_rho, 
               xlab ="Adult mosquito estimation (BG traps)" , 
               ylab = "Adult mosquito Simulations",
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

#--------------------------------------------------------------------------------------------#
################# CORRELATION BETWEEN RHOs ###################
Path_rho <- "~/MAD_MODEL/SUR_MODEL/Code/rho_sim.rds"
df_rho <- readRDS(Path_rho)
l = length(df_rho$rho)
df_rho$mov_rho <- mov_avg(l,dt,df_rho$rho)

Path_rho_ob <- "~/MAD_MODEL/SUR_MODEL/Code/rho_observed.rds"
df_rho_ob <- readRDS(Path_rho_ob)
l = length(df_rho_ob$rho)
df_rho_ob$mov_rho_1 <- mov_avg(l,dt,df_rho_ob$rho)
df_tot <- merge(df_rho, df_rho_ob, by = "date")

ggscatterstats(data = df_tot,
               x = mov_rho_1,
               y = mov_rho, 
               xlab ="Rho Observed" , 
               ylab = "Rho simulations") 
