rm(list = ls())
library(easypackages)
libraries("gdata", "ggExtra","ggplot2", "numbers","tidyverse","data.table","multiplex","reshape","viridis","stats","ggpubr","ggstatsplot","e1071","mlr3misc","deSolve", "gganimate") 

# STEPS to compute RHO(M).
# 1. Process the data to obtain A(t), registration temporal series, to do the simulations in C.
# 2. Take the output file Output_integration.dat from C.
# 3. Upload the reports file and the propensity probability.Check that all the time series goes from init_date 
# to end_date.
# 4. Compute Rho and plot it.

##### UPLOAD REGISTRATION FILE ###### STEP 1
# Data with registration for BCN taken from age distribution time series, filter by age 0.
Path_reg = "~/MAD_MODEL/MAD_MODEL/data/ages_days_bcn.csv"
registration = read.csv(Path_reg)
registration <- registration %>% filter( registration$age_days == 0)
registration$date <- as.Date(registration$date,'%Y-%m-%d') 
registration$age_days <- NULL
ref_date = min(registration$date)
# Remove the data after the update of the app.
init_date = min(registration$date)
end_date = max(registration$date)
registration$time = as.numeric(registration$date - as.Date(ref_date,"%Y-%m-%d") , units="days")
n_days = as.numeric(end_date - init_date) 
df_date <- data.frame(date = seq(init_date, end_date, "days"), time = seq(0,n_days,1))
df_test = data.frame(time = c(0:max(registration$time)), n =0)
# Añado los registros de los dias que no tengo para que tenga estos datos sin cortes.
registration <- merge(x = df_test, y = registration, by = "time", all.x = TRUE)
registration$N[is.na(registration$N)] <- 0
registration$n <- NULL
#registration$date <- NULL
downloads_bcn <- t(registration)

reg_plot <- merge(registration,df_date, by ="time")
ggplot(reg_plot) + 
  geom_line(aes(x = date.x, y= N), color = "darkblue")+
  ylab("Number of new users") +
  xlab("date") +
  scale_color_manual(values = c('#9E329F')) +
  theme_bw()+
  theme(text = element_text(size=19))

#---------------------------------------------------------------------------#
#### Data of participation from the simulation in C#####  STEP 2
Path <- "~/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_DETERMINISTIC/Output_Integration_bcn_2000_ages.data"
int_sol <- data.frame(t(read.table(Path, header=FALSE)))
# Adding age classes, to check the age distribution.
int_sum = rowSums(int_sol[,-1])
# Adding the first classes to check if there is a trend.
df_sum <- data.frame(time = int_sol[,1], Simulation = int_sum )
df_sum <- merge(df_sum,registration, by="time")
colnames(df_sum) <- c("time","Sum of all classes", "Registered")
df_plot <- reshape2::melt(df_sum, id.vars = c("time"))
df_plot <- merge(df_plot, df_date, by = "time")
# Number of participants age 1 vs registration.
df_aux <- int_sol[,c("X1","X2")]
colnames(df_aux) <- c("time","Age 1")
df_aux <- merge(df_aux, registration, by = "time")
df_plot <- reshape2::melt(df_aux, id.vars = c("time"))

#---------------------------------------------------------------------------#
########## PLOTS PARTICIPATION ###########
# Checking the sumation of the different age groups. 
mat_sim = as.matrix(int_sol)
int_sum_1 = mat_sim[,2:31]%*%prop_mat[1:30,4]
int_sum_2 = mat_sim[,31:649]%*%prop_mat[30:648,4]
int_sum_3 = mat_sim[,649:1386]%*%prop_mat[648:1385,4]
int_sum_4 = mat_sim[,951:1386]%*%prop_mat[950:1385,4]
int_sum = mat_sim[,2:1386]%*%prop_mat[1:1385,4]
df_sum <- data.frame(time = int_sol[,1], 
                     First_age_group = int_sum_1,
                     Second_age_group = int_sum_2,
                     Third_age_group = int_sum_3,
                     total  = int_sum,
                     registered = registration$N[1:length(int_sum)])
df_plot <- reshape2::melt(df_sum, id.vars = c("time"))
df_plot <- merge(df_plot, df_date, by = "time")

########## PLOTS PARTICIPATION ##############
# Plot to check which age group is increasing with the total sum and registered.
ggplot(df_plot,aes(date, value)) +
  geom_line(aes( colour = variable))  +
  scale_x_date(date_breaks = "10 month",
               date_labels = "%b %y")  +
  ylab("Counts") + ggtitle("Participation dynamics") +
  scale_color_manual(name = "",
                     labels = c("<30","[30:648]",
                                ">649", "Total","Registered, B(t)"),
                     values=c('#686868','#686868','#686868', '#0303cd','#686868'))+
  theme(text = element_text(size=16)) +
  theme_bw() 

df_plot <- df_plot %>% filter(variable != "total" )
ggplot(df_plot,aes(date, value)) +
  geom_line(aes( colour = variable))  +
  scale_x_date(date_breaks = "10 month",
               date_labels = "%b %y")  +
  ylab("Counts") + ggtitle("Participation dynamics") +
  scale_color_manual(name = "",
                     labels = c("<30","[30:648]",
                                ">649","Registered, B(t)"),
                     values=c('#c71585','#6495ed','#8fbc8f','#deb887'))+
  theme(text = element_text(size=16)) +
  theme_bw() 

df_plot <- df_plot %>% filter(variable == "registered" )
ggplot(df_plot,aes(date, value)) +
  geom_line(aes( colour = variable))  +
  scale_x_date(date_breaks = "10 month",
               date_labels = "%b %y")  +
  ylab("Counts") + ggtitle("Participation dynamics") +
  scale_color_manual(name = "",
                     labels = c("Registered, B(t)"),
                     values=c('#deb887'))+
  theme(text = element_text(size=16)) +
  theme_bw() 



# Plot with the registration dynamics.
df_sum <- data.frame(time = int_sol[,1], 
                     registered = registration$N[1:length(int_sum)])
df_plot <- reshape2::melt(df_sum, id.vars = c("time"))
df_plot <- merge(df_plot, df_date, by = "time")

ggplot(df_plot,aes(date, value)) +
  geom_line(aes( colour = variable))  +
  scale_x_date(date_breaks = "10 month",
               date_labels = "%b %y")  +
  ylab("Counts") + ggtitle("Participation dynamics") +
  scale_color_manual(name = "",
                     labels = c("Registered"),
                     values=c('#2F7D9F'))+
  theme(text = element_text(size=16))+
  theme_bw()

# Plot withn the registration dynamics, and the age groups.
df_sum <- data.frame(time = int_sol[,1], 
                     First_age_group = int_sum_1,
                     Second_age_group = int_sum_2,
                     Third_age_group = int_sum_3,
                     Fourth_age_group = int_sum_4,
                     #  total  = int_sum,
                     registered = registration$N[1:length(int_sum)])
df_plot <- reshape2::melt(df_sum, id.vars = c("time"))
df_plot <- merge(df_plot, df_date, by = "time")

ggplot(df_plot,aes(date, value)) +
  geom_line(aes( colour = variable))  +
  scale_x_date(date_breaks = "10 month",
               date_labels = "%b %y")  +
  ylab("Counts") + ggtitle("Participation dynamics") +
  scale_color_manual(name = "",
                     labels = c("[1:100]","[100:500]",
                                "[500:950]","[950:1385]","Registered"),
                     values=c('#FF00F6','#FF2C00','#00FF5E','#0092F6','#2F7D9F'))+
  theme(text = element_text(size=16))+
  theme_bw()

#---------------------------------------------------------------------------#
##### Moving Average #####
# Function to compute the moving average:
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

# Participation:
colnames(df_sum) <- c("time","Simulation", "Registered")
l = length(df_sum$time)
dt = 7
df_sum$mov = mov_avg(l,dt,df_sum$Simulation)
df_sum <- merge(df_sum, df_date, by = "time")
ggplot(df_sum) + 
  geom_line(aes(date, Simulation)) +
  scale_x_date(date_breaks = "6 month",
               date_labels = "%b %y", 
               limits = as.Date(c("2018-05-01","2020-12-01"))) +
  ylab("Moving average")  +
  ggtitle("Sumation of age classes")+
  theme_bw()
  
#---------------------------------------------------------------------------#
##### REPORTS UPLOAD ###### STEP 3
Path = "~/MAD_MODEL/MAD_MODEL/data/a000_mosquito_alert_spatio_temporal_data_D_mod_df.Rds"
# Path = paste(PC,Path, sep="")

reports = read_rds(Path) %>% filter(presence==TRUE)
reports$date= as.Date(reports$date,"%Y-%m-%d") 
reports$id = 1
reports$time = as.numeric(reports$date - as.Date(ref_date,"%Y-%m-%d") , units="days")
reports <- reports[,c("id","date", "time")]
min_date = min(reports$date)
max_date = max(reports$date)
# Group by date, this will create a table with each date and the number of reports
reports <- reports %>%  group_by(time) %>%  summarise(mean = mean(id), sum = sum(id), n = n())
reports <- reports[,c("time", "n")]
df_aux <- data.frame(time = c(min(reports$time):max(reports$time)), y = 1)
reports <- merge(df_aux, reports, by = "time", all.x = "TRUE")
reports <- reports[,c("time", "n")]
reports[is.na(reports)] <- 0
reports$n_smooth <- smooth(reports$n)
# Add the date to the df reports.
reports <- merge(reports,df_date, by ="time")

####### PLOT REPORTS #######

# Delete Outliers:
reports$n[reports$date > as.Date("2018-12-01" , "%Y-%m-%d") &
              reports$date < as.Date("2019-04-01" , "%Y-%m-%d")] <- 0
reports$n[reports$date > as.Date("2019-12-01" , "%Y-%m-%d") & 
              reports$date < as.Date("2020-04-01" , "%Y-%m-%d")] <- 0

ggplot(reports) + 
  geom_line(aes(date, n)) +
  #xlim(0,400) +
  xlab("Date") + 
  ylab("Number of reports") +
  scale_x_date(date_breaks = "6 month",
               date_labels = "%b %y")+
  scale_color_manual(values = c('#32329f')) +
  theme_bw()+
  theme(text = element_text(size=20))

#---------------------------------------------------------------------------#
###### PROPENSITY PROBABILITY UPLOAD#######
Path_prop <- "~/MAD_MODEL/MAD_MODEL/data/propensity_predictions.csv"
prop_mat  <- read.csv(Path_prop)

# File .dat of propensity probabilities.
prop_vec <- t(prop_mat[,4])
min_prop = min(prop_mat[,4])
max_prop = max(prop_mat[,4])
ggplot(prop_mat)+ 
  geom_line(aes(x =participation_time_days, y = reporting_prob))+
  xlab("Participants age (days)")+
  ylab("Probability")+
  scale_color_manual(values = c('#9E329F')) +
  theme_bw()+
  theme(text = element_text(size=21))

prop_mat$group <- 0
prop_mat$group[prop_mat$participation_time_days<=100]  <- "FIRST" 
prop_mat$group[prop_mat$participation_time_days>100 & prop_mat$participation_time_days<=500]  <- "SECOND" 
prop_mat$group[prop_mat$participation_time_days>500 & prop_mat$participation_time_days<=950]  <- "THIRD" 
prop_mat$group[prop_mat$participation_time_days>950]  <- "FOURTH" 
min_prop = min(prop_mat$participation_time_days)
max_prop = max(prop_mat$participation_time_days)

ggplot(prop_mat) + 
  geom_line(aes(participation_time_days,reporting_prob,color=group)) +
  xlab("Age of the participant (days)") + 
  ylab("Probability of reporting")+
  scale_x_continuous(breaks = round(seq(min_prop, max_prop, by = 100),1)) +
  scale_color_manual(values=c('#FF00F6','#FF2C00','#00FF5E','#0092F6'))+
  theme(text = element_text(size=16)) + theme_bw()


#---------------------------------------------------------------------------#
# Remove unused variables.

rm(df_plot);rm(df_aux);rm(df_sum)
rm(df_test); rm(rg_plot); rm(downloads_bcn)
rm(int_sum); rm(int_sum_1);rm(int_sum_2); rm(int_sum_3); rm(int_sum_4)

max_int = ncol(int_sol) - 1
max_prop = max(prop_mat[,1])
max_l = min(c(max_int, max_prop))
#---------------------------------------------------------------------------#
###### COMPUTE rho(M) ######## STEP 4
rho_t <- function(t, mat){
  mat <- as.matrix(mat)
  ind <- which(mat[,1] == t)
  den_sum <- mat[ind,2:1386]%*%prop_mat[1:1385,4]
  if(den_sum != 0){
    ind1 <- which(reports$time == t)
    if(length(ind1) == 0){
      rep <- 0
    }else{
      rep <- reports$n[ind1]
    }
    rho_m <- rep/den_sum
  }else{
    rho_m <- 0
  } 
  vec <- c(rho_m,den_sum)
  return(vec)
}

############# RHO COMPUTATION #############
# Simulated data:
t_init = max(c(min(int_sol$X1),min(reports$time)))
t_end = min(c(max(int_sol$X1),max(reports$time)))
vec = c(t_init:t_end)
rho = matrix(0,1,length(vec))
den = matrix(0,1,length(vec))
for(i in vec){
  output <- rho_t(i,int_sol)
  rho[i-t_init+1] = output[1]
  den[i-t_init+1] = output[2]
  print(paste0("i:",i))
}

df_rho <- data.frame(time = vec, rho = t(rho))
df_rho <- merge(df_rho, df_date, by = "time")

ggplot(df_rho) + 
  geom_line(aes(date, rho))  +
  scale_x_date(date_breaks = "3 month",date_labels = "%b %y") +
  ylab(expression(rho))+
  theme(text = element_text(size=14))+
  theme_bw()

filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/rho_sim.rds") 
saveRDS(df_rho, file = filename)

# Compute Rho with observed data:
Path = "~/MAD_MODEL/MAD_MODEL/data/Ob_data.rds"
ob_data <- readRDS(Path)
df_date_ob <- ob_data[,c(1,2)]
ob_data$date <- NULL
t_init = max(c(min(ob_data[,1]),min(reports$time)))
t_end = min(c(max(ob_data[,1]),max(reports$time)))
vec = c(t_init:t_end)
rho = matrix(0,1,length(vec))
den = matrix(0,1,length(vec))
for(i in vec){
  output <- rho_t(i,ob_data)
  rho[i-t_init+1] = output[1]
  den[i-t_init+1] = output[2]
  print(paste0("i:",i))
}

df_rho_ob <- data.frame(time = vec, rho = t(rho))
df_rho_ob <- merge(df_rho_ob, df_date_ob, by = "time")

ggplot(df_rho_ob) + 
  geom_line(aes(date, rho))  +
  scale_x_date(date_breaks = "3 month",date_labels = "%b %y") +
  ylab(expression(rho))+
  theme(text = element_text(size=14))+
  theme_bw()

filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/rho_observed.rds") 
saveRDS(df_rho_ob, file = filename)
#---------------------------------------------------------------------------# 

# Moving average to rho:
l = length(df_rho$rho)
dt = 7
mov = mov_avg(l,dt,df_rho$rho)
df_rho_mov <- data.frame(time = vec, rho = mov)
df_rho_mov <- merge(df_rho_mov, df_date, by = "time")
ggplot(df_rho_mov) + 
  geom_line(aes(date, rho))  +
  scale_x_date(date_breaks = "3 month",date_labels = "%b %y") +
  ylab(expression(rho))+
  theme(text = element_text(size=14))+
  theme_bw()

###### Upload data from bgtraps.######
Path = '~/MAD_MODEL/MAD_MODEL/data/bcn_bgcount_time_profile.csv'
# Path = paste(PC,Path, sep="")

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

## Maxima área de influencia a mano alzada: 500/(pi*300^2) .
ggplot(df_bg_mov) + 
  geom_line(aes(date, bg_counts))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab("Adults mosquitoes per square meter") +
  theme_bw()

# ###########SCATTER PLOTS########################
df_tot <- merge(df_rho_mov, df_bg_mov, by = "date")

# Scatter plot with marginal distributions.
ggscatterstats(data = df_tot,
               x = bg_counts,
               y = rho, 
               xlab ="Adult mosquito per m^2 (BG traps)" , 
               ylab =expression(rho),
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

# Scatter plot for the rho_M and BG traps estimates.
linear_mod <- lm(df_tot$rho~df_tot$bg_counts)
summary(linear_mod)

# Cleam memory:
rm(den)
rm(df_bg_mov)
rm(df_date)
rm(df_plot)
rm(df_rho)
rm(df_rho_mov)
rm(int_sol_fil)
rm(int_sum)
rm(int_sum_1)
rm(int_sum_2)
rm(int_sum_3)
rm(int_sum_4)
rm(mat_sim)
rm(ob_data)
