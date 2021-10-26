rm(list = ls())
library(easypackages)
libraries("scales","gdata", "ggplot2",
          "numbers","tidyverse","data.table",
          "multiplex","reshape","viridis",
          "stats","ggpubr","ggstatsplot",
          "e1071","mlr3misc","deSolve",
          "gganimate") 
# Different path for MAc and Ubuntu.
# PC = "/Users/celsaaraujobarja/Documents"
# PC = "/Users/celsaaraujobarja/Documents"
#PC = "/home/marta/Documentos"
#Path= "/home/marta/Documentos/PHD/2021/SUR_Model/Code/OUTPUT/df_rho.dat"
# Path_rel= "/PHD/2021/SUR_Model/Code/OUTPUT/df_rho.dat"
# Path = paste(PC,Path_rel, sep="")
Path = "~/MAD_MODEL/VECTOR_MODEL/data/df_rho.dat"
df_rho <- data.frame(t(read.table(Path, header=FALSE)))
colnames(df_rho) <- c("time", "rho", "date")
df_rho$date = as.Date(df_rho$date , "%Y-%m-%d")
df_rho$time <- NULL
df_rho$rho <- as.numeric(df_rho$rho)
ggplot(df_rho) + 
  geom_line(aes(x = date, y =rho)) +
  ggtitle("Encounter rate computed from Citizen Science model") +
  theme_bw()
# Temperatures from Barcelona:
Path_temp = "~/MAD_MODEL/VECTOR_MODEL/data/bcn_weather_daily.Rds"
# Path_temp = paste(PC,Path_temp, sep="")

temp <-read_rds(Path_temp)
temp$date = as.Date(temp$date , "%Y-%m-%d")
temp <- temp %>%  group_by(date) %>% summarise(mean_temp = mean(valor))

ggplot(temp) +
  geom_line(aes(date, mean_temp))+
  ggtitle("Mean temperature Barcelona") + 
  xlab("Mean temperature")+
  theme_bw()

#########GONOTROPHIC CYCLE#########
temp_vec = c(20,25,30,35)
dur_gono = c(6.7,4,2.9,4.7)
temp_vec2 <- temp_vec^2
quadratic_mod <- lm(dur_gono ~ temp_vec + temp_vec2)
temp_val <- seq(0,40,1)
predict_gono <- predict(quadratic_mod, list(temp_vec=temp_val, temp_vec2<- temp_val^2))
df_gono_data <- data.frame(temp_vec, dur_gono)
df_predict <- data.frame(temp_val, predict_gono)

ggplot(df_predict, aes(temp_val, predict_gono)) +
  geom_line()+
  geom_point(data = df_gono_data, aes(x = temp_vec, y =dur_gono))+
  theme_bw()

int = quadratic_mod$coefficients[1]
beta = quadratic_mod$coefficients[2]
beta2 = quadratic_mod$coefficients[3]
int
beta
beta2
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
# df_dl_opt_vec <- data.frame(temp = vec, deltaA = unlist(lapply(vec,dL_opt)))

ggplot(df_gonot_vec) + geom_line(aes(temp,gonot)) + 
  ggtitle("Inverse of the Gonotrophic cycle")+
  theme_bw()

ggplot(df_dL_vec) + geom_line(aes(temp,dL)) +
  ggtitle("Larva development rate")+
  theme_bw()

ggplot(df_deltaL_vec) + geom_line(aes(temp,deltaL)) +
  ggtitle("Larva mortality rate")+
  theme_bw()

ggplot(df_deltaA_vec) + geom_line(aes(temp,deltaA)) +
  ggtitle("Mosquito adult mortality rate")+
  theme_bw()

ggplot(df_dl_opt_vec) + geom_line(aes(temp,deltaA)) +
  ggtitle("Mosquito adult mortality rate")+
  theme_bw()


# Compute the values of the functions/forcings with temp.
# Compute the minimum date of the rho:
min_rho_date <- min(df_rho$date)
min_temp_date <- min(temp$date)
min_date <-max(min_rho_date, min_temp_date)

max_rho_date <- max(df_rho$date)
max_temp_date <- max(temp$date)
max_date <-min(max_rho_date, max_temp_date)
# DFs with the date and value of the parameter at that time.
temp <- temp  %>% filter( temp$date >= min_date & temp$date <= max_date)
df_date <- data.frame(date = temp$date)
df_date$time = as.numeric(df_date$date - as.Date(min_date,"%Y-%m-%d") , units="days") 

gono = unlist(lapply(temp$mean_temp,gonot))
df_gonot_out <- data.frame(date = temp$date, gono)

df_gonot_out$time = as.numeric(df_gonot_out$date - as.Date(min_date,"%Y-%m-%d") , units="days") 
df_gonot_out <- df_gonot_out %>% filter( df_gonot_out$time >= 0)

df_dL_out <- data.frame(date = temp$date, dL = unlist(lapply(temp$mean_temp,d_L)))
df_dL_out$time = as.numeric(df_dL_out$date - as.Date(min_date,"%Y-%m-%d") , units="days") 
df_dL_out <- df_dL_out %>% filter( df_dL_out$time >= 0)

df_deltaL_out <- data.frame(date = temp$date, deltaL = unlist(lapply(temp$mean_temp,delta_L)))
df_deltaL_out$time = as.numeric(df_deltaL_out$date - as.Date(min_date,"%Y-%m-%d") , units="days") 
df_deltaL_out <- df_deltaL_out %>% filter( df_deltaL_out$time >= 0)

df_deltaA_out <- data.frame(date = temp$date, deltaA = unlist(lapply(temp$mean_temp,delta_A)))
df_deltaA_out$time = as.numeric(df_deltaA_out$date - as.Date(min_date,"%Y-%m-%d") , units="days") 
df_deltaA_out <- df_deltaA_out %>% filter( df_deltaA_out$time >= 0)

df_rho$time = as.numeric(df_rho$date - as.Date(min_date,"%Y-%m-%d") , units="days") 
df_rho <- df_rho %>% filter( df_rho$time >= 0)

ggplot(df_gonot_out) +
  geom_line(aes(date,gono)) +
  ggtitle("Gonotrophic cycle")+
  theme_bw()

ggplot(df_dL_out) +
  geom_line(aes(date,dL)) +
  ggtitle("Larva development rate")+
  theme_bw()

ggplot(df_deltaL_out) +
  geom_line(aes(date,deltaL)) +
  ggtitle("Larva mortality rate")+
  theme_bw()

ggplot(df_deltaA_out) +
  geom_line(aes(date,deltaA)) +
  ggtitle("Adult mosquito mortality rate")+
  theme_bw()

df_gonot_out$date <- NULL
df_gonot_out <- df_gonot_out[,c(2,1)]
df_dL_out$date <- NULL
df_dL_out <- df_dL_out[,c(2,1)]
df_deltaL_out$date <- NULL
df_deltaL_out <- df_deltaL_out[,c(2,1)]
df_deltaA_out$date <- NULL
df_deltaA_out <- df_deltaA_out[,c(2,1)]
df_rho$date <- NULL
df_rho <- df_rho[,c(2,1)]


###############   ODE INTEGRATION   ##################
# require(deSolve)
# library.dynam.unload("deSolve", libpath=paste(.libPaths()[1], "//deSolve", sep=""))
# library.dynam("deSolve", package="deSolve", lib.loc=.libPaths()[1])
# OJOOOOO!!! Cuando cambias de PC borrar .o y .so.
Path = "~/MAD_MODEL/VECTOR_MODEL/Code/INTEGRATION/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model.c")
dyn.load("model.so")

f = 200
K = 250000
H = 1600000
# We create a vector with the constant parameters.
parms = c(f,K,H)
# We set the initial conditions to cero.
Y <- c(y1 = 100.0, y2 = 0.0, y3 = 0.0)
# List with the data frames of the forcings, sort as the c code.
forcs_mat <- list(data.matrix(df_gonot_out),
                  data.matrix(df_dL_out),
                  data.matrix(df_deltaL_out),
                  data.matrix(df_rho),
                  data.matrix(df_deltaA_out))
min_t <- min(df_rho$time)
max_t <- max(df_rho$time)
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model",
           initfunc = "initmod", nout = 1,
           outnames = "Sum", initforc = "forcc",
           forcings = forcs_mat, fcontrol = list(method = "constant")) 

ode <- data.frame(out) 
ode_df <- merge(ode, df_date, by ="time")
ode_df$Sum <- NULL
head(ode_df)
colnames(ode_df) <- c("Time", "L", "Ah","A", "date" )
df_L_A <- ode_df[,c(5,2,3)]
df_Ah <- ode_df[,c(5,4)]
df_L <- ode_df[,c(5,2)]
df_plot_1 <- reshape2::melt(df_L_A, id.vars = c("date"))
df_plot_L <- reshape2::melt(df_L, id.vars = c("date"))
df_plot <- reshape2::melt(df_Ah, id.vars = c("date"))

scientific_10 <- function(x) {
  parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x)))
}

ggplot(df_plot_1,aes(date, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Counts") +
  ggtitle("Vector dynamics")+
  scale_color_manual(name = "",
                     labels = c("Larva", "Adult mosquito"),
                     values=c('#FF00F6','#FF2C00'))+
  theme_bw() + scale_y_continuous(labels = scientific_10)+
  theme(text = element_text(size=18))

ggplot(df_Ah)  +
  geom_line(aes(date, A), color = "dark green") +
  ylab("Counts") +
  ggtitle("Handling mosquitoes dynamics") +
  scale_color_manual(values=c('#FF00F6')) +
  theme_bw() + scale_y_continuous(labels = scientific_10)+
  theme(text = element_text(size=18))

ggplot(df_L)  +
  geom_line(aes(date, L), color = "dark green") +
  ylab("Counts") +
  ggtitle("Larvae dynamics") +
  scale_color_manual(values=c('#FF00F6')) +
  theme_bw()+
  theme(text = element_text(size=16))

###### Equilibrium points ######

eq_point<- function(a,f,dL,deltaA,deltaL,chi,H,K){
  A_h <- K*((chi*H*dL/((a+deltaA)*(chi*H+deltaA)))-(1/(a*f))*(dL+deltaL))
  L <- (a*f*A_h)/(((a*f*A_h)/K)+(dL+deltaL))
  A <- (dL*K/(chi*H+deltaA))*L
  return(c(L,A_h,A))
}

gono_vec <- df_gonot_out$gono
dL_vec <- df_dL_out$dL
deltaL_vec <- df_deltaL_out$deltaL
deltaA_vec <- df_deltaA_out$deltaA

chi <- 0.0114512
  
eq_out <- mapply(eq_point, gono_vec,f,dL_vec,deltaA_vec,deltaL_vec,chi,H,K)

eq_df <- data.frame(time = df_date$date, t(eq_out))
colnames(eq_df) <- c("time","L","A_h","A")
eq_df <- reshape2::melt(eq_df, id.vars = c("time"))

ggplot(eq_df,aes(time, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Number of individuals") +
  ggtitle("Equilibrium point")+
  scale_color_manual(name = "",
                     labels = c("Larva", "Handling adult mosquito", "Adult mosquito"),
                     values=c('#FF00F6','#FF2C00', '#FF2C00')) +
  theme_bw()
