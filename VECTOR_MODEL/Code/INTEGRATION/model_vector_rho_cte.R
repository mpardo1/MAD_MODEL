rm(list = ls())
library(easypackages)
libraries("scales","gdata", "ggplot2",
          "numbers","tidyverse","data.table",
          "multiplex","reshape","viridis",
          "stats","ggpubr","ggstatsplot",
          "e1071","mlr3misc","deSolve",
          "gganimate") 
# Rho data:
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



###############   ODE INTEGRATION   ##################

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
