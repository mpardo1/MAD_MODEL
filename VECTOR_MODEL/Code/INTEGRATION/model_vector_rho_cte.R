rm(list = ls())
library(easypackages)
libraries("ggplot2","tidyverse","ggstatsplot","deSolve")

# Rho data:
Path = "~/MAD_MODEL/VECTOR_MODEL/data/df_rho.dat"
df_rho <- data.frame(t(read.table(Path, header=FALSE)))
colnames(df_rho) <- c("time", "rho", "date")
df_rho$date = as.Date(df_rho$date , "%Y-%m-%d")
df_rho$rho <- as.numeric(df_rho$rho)

min_date <- min(df_rho$date)
max_date <- max(df_rho$date)
df_rho$time <- as.numeric(df_rho$date - as.Date(min_date,"%Y-%m-%d") , units="days")
df_rho$rho[df_rho$date > as.Date("2018-12-01" , "%Y-%m-%d") &
             df_rho$date < as.Date("2019-04-01" , "%Y-%m-%d")] <- 0
df_rho$rho[df_rho$date > as.Date("2019-12-01" , "%Y-%m-%d") & 
             df_rho$date < as.Date("2020-04-01" , "%Y-%m-%d")] <- 0

ggplot(df_rho) + 
  geom_line(aes(x = date, y =rho)) +
  ggtitle("Encounter rate computed from Citizen Science model") +
  theme_bw() +
  theme(text = element_text(size=14))

df_date <- df_rho[,c(1,3)]
###############   ODE INTEGRATION   ##################

Path = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_vec_cte.c")
dyn.load("model_vec_cte.so")

fec = 100
K = 250000
Hum = 13554
omega_t = 0.2
delta_L = 0.2
delta_A = 0.3
d_L = 0.8
a = 0.01

parms = c(fecun = fec, Ka = K, Hu = Hum, del_L = delta_L, del_A = delta_A, dev_L = d_L, gon = a)
# We set the initial conditions to zero.
Y <- c(y1 = 100, y2 = 0, y3 = 0)
min_t <- 1
max_t <- 100
times <- seq(min_t,max_t, 1)
# List with the data frames of the forcings, sort as the c code.
# df_rho$rho[df_rho$rho == 0 ] <- 0.00000001
forcs_mat <-data.matrix(df_rho[,1:2])

min_t <- min(df_rho$time)
max_t <- max(df_rho$time)
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model_vec_cte", method = "ode45",
           initfunc = "initmod", nout = 1,
           outnames = "Sum", initforc = "forcc",
           forcings = forcs_mat, fcontrol = list(method = "constant")) 

# Event :
# This is where we define your event function
# Add this directly above your call to ode()
posfun <- function(t, y, parms){
  with(as.list(y), {
    y[which(y<0)] <- 0  
    return(y)
  })
}


out <- ode(Y, times=times, func = "derivs",
           parms = parms, dllname = "model_vec_cte", method = "ode45",
           initfunc = "initmod", nout = 1,
           outnames = "Sum", initforc = "forcc",
           forcings = forcs_mat, fcontrol = list(method = "constant"),
           events=list(func = posfun, time = c(0:max_t)))

ode <- data.frame(out) 

ode_df <- merge(ode, df_date, by ="time")
ode_df$Sum <- NULL
head(ode_df)
colnames(ode_df) <- c("Time", "L", "Ah","A", "date" )
df_L_A <- ode_df[,c(5,2,3)]
df_Ah <- ode_df[,c(5,3)]
df_L <- ode_df[,c(5,2)]
df_A <- ode_df[,c(5,4)]
df_plot_1 <- reshape2::melt(df_L_A, id.vars = c("date"))
df_plot_L <- reshape2::melt(df_L, id.vars = c("date"))
df_plot_Ah <- reshape2::melt(df_A, id.vars = c("date"))
df_plot <- reshape2::melt(df_Ah, id.vars = c("date"))
df_A$
saveRDS(df_A, file = "~/MAD_MODEL/SUR_MODEL/Code/adults.rds")

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

ggplot(df_L)  +
  geom_line(aes(date, L), color = "dark green") +
  ylab("Counts") +
  ggtitle("Larvae dynamics") +
  scale_color_manual(values=c('#FF00F6')) +
  theme_bw()+
  theme(text = element_text(size=16))

ggplot(df_A)  +
  geom_line(aes(date, A), color = "dark green") +
  ylab("Counts") +
  ggtitle("Adult mosquito dynamics") +
  scale_color_manual(values=c('#FF2C00')) +
  theme_bw()+
  theme(text = element_text(size=16))


ggplot(df_Ah)  +
  geom_line(aes(date, Ah), color = "dark green") +
  ylab("Counts") +
  ggtitle("Handling Adult mosquito dynamics") +
  theme_bw()+
  theme(text = element_text(size=16)) +
  scale_color_manual(values=c('#f0ff33')) 
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
