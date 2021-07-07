rm(list = ls())
library(easypackages)
libraries("gdata", "ggExtra","ggplot2", "numbers","tidyverse","data.table","multiplex","reshape","viridis","stats","ggpubr","ggstatsplot","e1071","mlr3misc","deSolve", "gganimate") 

# STEPS to compute RHO(M).
# 1. Process the data to obtain A(t), registration temporal series, to do the simulations in C.
# 2. Take the output file Output_integration.dat from C.
# 3. Upload the reports file and the propensity probability.Check that all the time series goes from init_date 
# to end_date.
# 4. Compute Rho and plot it.


# Different path for MAc and Ubuntu.
# PC = "/Users/celsaaraujobarja/Documents"
#PC = "/home/marta/Documentos"
##### UPLOAD REGISTRATION FILE ###### STEP 1
# Data with registration for BCN taken from age distribution time series, filter by age 0.
Path_reg = "~/MAD_MODEL/MAD_MODEL/data/ages_days_bcn.csv"
# Path_reg = paste(PC,Path_reg, sep="")

registration = read.csv(Path_reg)
registration <- registration %>% filter( registration$age_days == 0)
registration$date <- as.Date(registration$date,'%Y-%m-%d') 
registration$age_days <- NULL
ref_date = min(registration$date)
# Remove the data after the update of the app.
init_date = min(registration$date)
end_date = max(registration$date)
# registration <- registration %>% filter(date < as.Date(end_date,"%Y-%m-%d") )
# registration <- registration %>% filter(date >= as.Date(init_date,"%Y-%m-%d") )
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
# +
#   ggtitle("Registration time series Barcelona")
# Save FILE with DOWNLOADS DATA, A(t) en el modelo. This file will be save with registration.dat name.
# write.dat(downloads_bcn, "/home/marta/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_DETERMINISTIC/")

#### Data of participation from the simulation in C#####  STEP 2
PC_1 = "/Users/celsaaraujobarja"
#PC_1 = "/home/marta"
Path = "~/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_DETERMINISTIC/Output_Integration_bcn_2000_ages.data"
# Path = paste(PC_1,Path, sep="")

# Mac PATH
#Path = "/Users/celsaaraujobarja/Documents/PhD/Output_Int/Output_Integration.dat"
int_sol = data.frame(t(read.table(Path, header=FALSE)))
int_sol_fil = int_sol[,c(1:1000)]
int_sol_fil <- int_sol
int_sol_fil <- reshape2::melt(int_sol_fil, id.vars = c("X1"))
int_sol_fil$age_class <-as.numeric(sub('.','',int_sol_fil$variable))
ggplot(int_sol_fil) +
  geom_point(aes(age_class,value)) +
  transition_time(X1)+
  theme_bw()

# write.dat(int_sol[1:20,1:20], "/home/marta/Documentos/PHD/2021/SUR_Model/Code/OUTPUT/")
# write.dat(int_sol, "/home/marta/Documentos/PHD/2021/SUR_Model/Code/OUTPUT/")

# Adding age classes, to check the age distribution.
int_sum = rowSums(int_sol[,-1])
# Adding the first classes to check if there is a trend.
# int_sum = rowSums(int_sol[1:2300,2:30])
df_sum <- data.frame(time = int_sol[,1], Simulation = int_sum )
df_sum <- merge(df_sum,registration, by="time")
colnames(df_sum) <- c("time","Sum of all classes", "Registered")
df_plot <- reshape2::melt(df_sum, id.vars = c("time"))
df_plot <- merge(df_plot, df_date, by = "time")


# Plot number of participants at each time vs registration.
ggplot(df_plot,aes(date, value)) +
scale_x_date(date_breaks = "6 month",
  date_labels = "%b %y", 
  limits = as.Date(c("2018-05-01","2020-12-01"))) +
geom_line(aes( colour = variable)) +
ylab("Counts") +
scale_color_manual(values = c('#9E329F','#1642FE'))+
  ggtitle("Participation dynamics")+
  theme_bw()

ggplot(df_plot,aes(date, value)) +
  scale_x_date(date_breaks = "6 month",
               date_labels = "%b %y") +
  geom_line(aes( colour = variable)) +
  scale_color_manual(values = c('#9E329F','#1642FE')) +
  ggtitle("Participation dynamics")+
  theme_bw()
# Number of participants age 1 vs registration.
df_aux <- int_sol[,c("X1","X2")]
colnames(df_aux) <- c("time","Age 1")
df_aux <- merge(df_aux, registration, by = "time")
df_plot <- reshape2::melt(df_aux, id.vars = c("time"))

ggplot(df_plot,aes(time, value)) + 
  geom_line(aes( colour = variable))  +
  ylab("Number of participants") +
#  xlim(0,10) +
  scale_color_manual(values=c('#1642FE','#9E329F')) +
  scale_linetype_manual(values = c("solid", "dotted"))+
  theme_bw()


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
  

##### REPORTS UPLOAD ###### STEP 3
Path = "~/MAD_MODEL/MAD_MODEL/data/a000_mosquito_alert_spatio_temporal_data_D_mod_df.Rds"
# Path = paste(PC,Path, sep="")

reports = read_rds(Path) %>% filter(presence==TRUE)
reports$date= as.Date(reports$date,"%Y-%m-%d") 
reports$id = 1
reports$time = as.numeric(reports$date - as.Date(ref_date,"%Y-%m-%d") , units="days")
# Filter the data to have the same dates for the registration A(t) and the reporting.
#reports <- reports %>% filter(date < as.Date(end_date,"%Y-%m-%d") )
#reports <- reports %>% filter(date >= as.Date(init_date,"%Y-%m-%d") )
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
report_fil = reports
ggplot(report_fil) + 
  geom_line(aes(date, n)) +
  #xlim(0,400) +
  xlab("Date") + 
  ylab("Number of reports") +
  scale_x_date(date_breaks = "6 month",
               date_labels = "%b %y")+
  scale_color_manual(values = c('#32329f')) +
  theme_bw()+
  theme(text = element_text(size=20))
# +
#   ggtitle(expression("Tiger mosquito reports to MA in Barcelona"))



# # Total Mosquito:
# l = length(reports$date)
# reports$mov_tot = 0.6
# dt = 7
# reports$mov_tot = mov_avg(l,dt,reports$Total_Mosquitos)
# 
# # Register:
# l = length(reg_group_sort$n)
# reg_group_sort$mov_n = 0.6
# dt = 7
# reg_group_sort$mov_n = mov_avg(l,dt,reg_group_sort$n)
# 
# # Plot Moving average of reporting.
# ggplot(reports) + 
#   geom_line(aes(date, mov)) +
#   #xlim(0,400) +
#   xlab("Date") + 
#   ylab("Moving Average")
# # Plot Moving average of reporting filter by year.
# ggplot(reports[1:365,]) + 
#   geom_line(aes(Date, mov)) +
#   #xlim(0,400) +
#   xlab("Date") + 
#   ylab("Moving Average")

###### PROPENSITY PROBABILITY UPLOAD#######
#Path_prop = "/Users/celsaaraujobarja/Documents/PhD/Output_Int/propensity_predictions.csv"
Path_prop= "~/MAD_MODEL/MAD_MODEL/data/propensity_predictions.csv"
# Path_prop = paste(PC,Path_prop, sep="")

prop_mat = read.csv(Path_prop)

# File .dat of propensity probabilities.
prop_vec <- t(prop_mat[,4])
# write.dat(prop_vec, "/home/marta/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_DETERMINISTIC")

min_prop = min(prop_mat[,4])
max_prop = max(prop_mat[,4])
ggplot(prop_mat)+ 
  geom_line(aes(x =participation_time_days, y = reporting_prob))+
  # scale_x_continuous(breaks = round(seq(min_prop, max_prop, by = 250),1))+
  xlab("Participants age (days)")+
  ylab("Probability")+
  scale_color_manual(values = c('#9E329F')) +
  theme_bw()+
  theme(text = element_text(size=21))
# +
#   ggtitle("Propensity probability")
# Plot propensity probability.
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

# write.dat(prop_mat[1:20,c(1,4)], "/home/marta/Documentos/PHD/2021/SUR_Model/Code/OUTPUT/")
# Remove unused variables.

rm(reg_group)
rm(df_aux)
rm(df_plot)
rm(df_sum)
rm(df_test)
rm(report_fil)

max_int = ncol(int_sol) - 1
max_prop = max(prop_mat[,1])
max_l = min(c(max_int, max_prop))
#max_l = 30
###### COMPUTE rho(M) ######## STEP 4
rho_t <-function(t, p, mat){
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

t_init = max(c(min(int_sol$X1),min(reports$time)))
t_end = min(c(max(int_sol$X1),max(reports$time)))
vec = c(t_init:t_end)
rho = matrix(0,1,length(vec))
den = matrix(0,1,length(vec))
p = 1
for(i in vec){
  output <- rho_t(i,p,int_sol)
  rho[i-t_init+1] = output[1]
  den[i-t_init+1] = output[2]
  print(paste0("i:",i))
}

df_rho <- data.frame(time = vec, rho = t(rho))
df_rho <- merge(df_rho, df_date, by = "time")

# Create data frame input data C:
df_input_date <- merge(df_rho, registration, by = "date")
min_date <- min(df_input_date$date)
max_date <- max(df_input_date$date)
df_input_date$time = as.numeric(df_input_date$date - as.Date(min_date,"%Y-%m-%d") , units="days")
diff = as.numeric(max_date - min_date, units="days")
vec_temp <- seq(0,diff,1)
df_aux <- data.frame(time = vec_temp, a = 0)
df_input_date <- merge(df_input_date,df_aux, by ="time", all.y = "TRUE")
df_input_date$rho[is.na(df_input_date$rho)] <- 0
df_input_date$N[is.na(df_input_date$N)] <- 0
df_input_date$date <- NULL
df_input_date$a.x <- NULL
df_input_date$a.y <- NULL
df_input_date <- df_input_date[,c("time","N","rho")]
# write.dat(t(df_input_date), "/home/marta/PROJECT_MOSQUITO_ALERT/MODEL_CALCULATIONS/TEMPORAL_EVOLUTION_DETERMINISTIC/")

df_rho$rho.X1 <- NULL
colnames(df_rho) <- c("time","rho","date")
df_rho <- t(df_rho)
# write.dat(df_rho, "/home/marta/Documentos/PHD/2021/SUR_Model/Code/OUTPUT/")

df_rho <- data.frame(time = vec, rho = t(rho), den = t(den))
df_merg <- reshape2::melt(df_rho, id.vars = "time")

# Plot rho and denominator.
ggplot(df_merg,aes(time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F','#1642FE'))+
  theme_bw() +theme(text = element_text(size=16))+
  ggtitle("Registration time series Barcelona")

ggplot(df_rho) + 
  geom_line(aes(time, den))  +
  scale_color_manual(values=c('F5329F','#1642FE'))+
  theme_bw() +theme(text = element_text(size=16))+
  ggtitle("Registration time series Barcelona")


df_rho$smooth_rho <- smooth(df_rho$rho)
df_tot <- merge(df_rho, reports, by="time")
# Plot rho(M,t)

ggplot(df_tot) + 
  geom_line(aes(date, smooth_rho))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab(expression(rho))+
  theme_bw()

# Plot rho(M,t) smooth.
ggplot(df_tot) + 
  geom_line(aes(date,rho))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab(expression(rho))+
  theme_bw()+theme(text = element_text(size=18))
# +
#   ggtitle("Encounter per human capita simulation")

ggplot(df_tot) + 
  geom_line(aes(date, n_smooth)) +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y")+
  ylab("Number of reports")+
  theme_bw()+theme(text = element_text(size=18))+
  ggtitle("Registration time series Barcelona")


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

l = length(bg_traps$prop)
dt = 7
mov = mov_avg(l,dt,bg_traps$prop)
df_bg_mov <- data.frame(date = bg_traps$date, bg_counts = mov)

## Maxima área de influencia a mano alzada: 500/(pi*300^2) .
ggplot(df_bg_mov) + 
  geom_line(aes(date, bg_counts))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab("Adults mosquitoes per square meter") +
  theme_bw()

l = length(bg_traps$prop)
dt = 7
mov = mov_avg(l,dt,bg_traps$prop)
df_bg_mov <- data.frame(date = bg_traps$date, bg_counts = mov)

df_tot <- merge(df_rho_mov, df_bg_mov, by = "date")
df_tot$time <- NULL
df_plot <- reshape2::melt(df_tot, id.vars = c("date"))
# Plot moving average rho and adult mosquitoes per human from BG
ggplot(df_plot,aes(date,value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F','#1642FE'))


###########SCATTER PLOTS########################
# Plot of the correlation between the BG estimation and the rho(M) all MOVING AVERAGE.
l = length(bg_traps$prop)
dt = 7
mov = mov_avg(l,dt,bg_traps$prop)
df_bg_mov <- data.frame(date = bg_traps$date, mov_bg_count = mov)

l = length(df_rho$rho)
dt = 7
mov = mov_avg(l,dt,df_rho$rho)
df_rho_mov <- data.frame(time = vec, rho = mov)
df_rho_mov <- merge(df_rho_mov, df_date, by = "time")
df_tot <- merge(df_rho_mov, df_bg_mov, by = "date")

# Scatter plot with marginal distributions.
ggscatterstats(data = df_tot,
               x = mov_bg_count,
               y = rho, 
               xlab ="Adult mosquito per m^2 (BG traps)" , 
               ylab =expression(rho),
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 


ggscatterstats(
  data = df_tot,
  x = mov_bg_count,
  y = rho,
  # making further customization with `ggplot2` functions
  ggplot.component = list(ggplot2::geom_rug(sides = "b"))
)
# Scatter splot with Pearson correlation coeficient.
ggscatter(df_tot, x = "mov_bg_count", y = "rho", 
          add = "reg.line", conf.int = TRUE, 
          color = "black", shape = 20, size = 1 ,
          add.params = list(color = "blue", fill = "darkblue"),
          cor.coef = TRUE, cor.method = "pearson",
          xlab ="Adult mosquito per m^2 (BG traps)" , 
          ylab =paste0(expression(rho),"(M,t)") ) + 
  ggtitle("Moving average weekly")+
  theme(text = element_text(size=16))

# Scatter plot for the rho_M and BG traps estimates.
# df_rho <- merge(df_rho, df_date, by = "time")
# df_tot <- merge(bg_traps, df_rho, by = "date")
ggscatter(df_tot, x = "rho", y = "mov_bg_count", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "pearson",
          xlab = paste0(expression(rho),"(M,t)"), ylab = "Adult mosquito per human (BG traps)")

scatter.smooth( x = df_tot$mov_bg_count , y = df_tot$rho) + 
  xlab("Counts BG") + 
  ylab(expression(rho))

linear_mod <- lm(df_tot$rho~df_tot$mov_bg_count)
summary(linear_mod)
#-----------------------------------------------------------------------------------------------
#########RHO OBSERVED#############
Path_reg = "~/MAD_MODEL/MAD_MODEL/data/ages_days_bcn.csv"
# Path_reg = paste(PC,Path_reg, sep="")

ages =  read.csv(Path_reg)
# Convert to data type date the registration time.
ages$date = as.Date(ages$date,"%Y-%m-%d") 
# Remove the data after the update of the app.
ages <- ages %>% filter(date < as.Date("2020-10-01","%Y-%m-%d") )
ages$time = as.numeric(ages$date - as.Date(ref_date,"%Y-%m-%d") , units="days") + 1
time <- ages$time
ages$age_days = ages$age_days + 1
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
mat_ages = t(mat_ages)
dim = max(ages$time)
mat_ages[1,] = c(0:dim)
mat_ages = t(mat_ages)
# Computation of Rho with observed data of participation:
t_init = max(c(min(mat_ages[,1]),min(reports$time)))
t_end = min(c(max(mat_ages[,1]),max(reports$time)))
vec = c(t_init:t_end)
rho = matrix(0,1,length(reports$time)-1)
den = matrix(0,1,length(reports$time)-1)
p = 1
max_int = ncol(mat_ages) - 1
max_prop = max(prop_mat[,1])
max_l = min(c(max_int, max_prop))
for(i in vec){
  rho[i-t_init+1] = rho_t(i,p,mat_ages)[1]
  den[i-t_init+1] = rho_t(i,p,mat_ages)[2]
  print(paste0("i-t_init:",i-t_init))
  #print(paste0("i:",i))
  print(paste0("den_sum[i]:",rho_t(i,p,mat_ages)[2]))
}

df_rho = data.frame(time = vec, rho = rho, den = den)
df_merg <- reshape2::melt(df_rho, id.vars = "time")

# Plot rho and denominator.
ggplot(df_merg,aes(time, value)) + 
  geom_line(aes( colour = variable))  +
  scale_color_manual(values=c('#9E329F','#1642FE'))


df_rho$smooth_rho <- smooth(df_rho$rho)
df_tot <- merge(df_rho, reports, by="time")
# Plot rho(M,t)

ggplot(df_tot) + 
  geom_line(aes(date, smooth_rho))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab("rho(A,t)")

# Plot rho(M,t) smooth.
ggplot(df_tot) + 
  geom_line(aes(date,rho))  +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y") +
  ylab("rho(A,t)")

ggplot(df_tot) + 
  geom_line(aes(date, n_smooth)) +
  scale_x_date(date_breaks = "6 month",date_labels = "%b %y")+
  ylab("Number of reports")


# Sum of the participants.
int_sum = rowSums(mat_ages[1:2300,2:30])
df_sum <- data.frame(time = mat_ages[1:2300,1], Simulation = int_sum )
df_sum <- merge(df_sum,registration, by="time")
colnames(df_sum) <- c("time","Simulation", "Registered")
df_plot <- reshape2::melt(df_sum, id.vars = c("time"))
df_plot <- merge(df_plot, df_date, by = "time")


# Plot number of participants at each time vs registration.
ggplot(df_plot,aes(date, value)) +
  scale_x_date(date_breaks = "6 month",
               date_labels = "%b %y", 
               limits = as.Date(c("2018-05-01","2020-12-01"))) +
  geom_line(aes( colour = variable)) +
  ylim(0,200) + 
  ylab("Counts") +
  scale_color_manual(values = c('#9E329F','#1642FE'))
# Cleam memory:
rm(int_sol)
rm(mat_ages)
rm(mat_fil)
rm(mat_t)
rm(report_fil)
rm(init_sol_mat)
rm(int_sol_t)
rm(mat_fil_t)


