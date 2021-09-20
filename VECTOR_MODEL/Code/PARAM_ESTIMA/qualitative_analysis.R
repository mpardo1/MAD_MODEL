rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("coda")

######-------------------------------FUNCTIONS--------------------------------######
# Equilibrium points:
eq_point <- function(delta_L, delta_A, d_L, a, fec, K, Hum, omega_t){
  Ah_eq <-  K*(((omega_t*Hum*d_L)/((a+delta_A)*(omega_t*Hum +delta_A)))-((1/(a*fec))*(d_L+delta_L)))
  L_eq <- ((a*fec*Ah_eq)/((a*fec/K)*Ah_eq+(d_L+delta_L)))
  A_eq <- ((d_L/(omega_t*Hum+delta_A))*((a*fec*Ah_eq)/(((a*fec/K)*Ah_eq)+(d_L+delta_L))))
  eq <- c(L_eq, A_eq, Ah_eq)
  return(eq)
}


ode_value <- function(L, A, Ah){
  val <- c()
  val[1] = a*fec*Ah*(1-(L/K))-(d_L+delta_L)*L
  val[2] = d_L*L - (omega_t*Hum + delta_A)*A
  val[3] = omega_t*Hum*A - (a + delta_A)*Ah
  return(val)
}

# Feasability condition:
feasability_cond <- function(del_L, del_A, dev_L, gon, fec, Kar, Hum, omega){
  fes <- 0
  val <- (gon*fec*omega*Hum*dev_L)/((gon + del_A)*(omega*Hum + del_A)*(dev_L + del_L))
  if(val > 1){
    fes <- 1
  }
  return(fes)
}

omega_stab <- function(omega){
  fes <- 0
  val <- delta_A/(((a*fec*Hum*d_L)/((a+delta_A)*(d_L+delta_L)))-Hum)
  if(omega > val ){
    fes <- 1
  }
  return(fes)
}

coeff_omega <- function(omega){
  val <- Hum*omega/(Hum*omega + delta_A)
  return(val)
}

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
###############   ODE INTEGRATION   ##################
#------------------------------------------------------------------------------#
# ODE system in R:
vect <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    
    dL <-  gon*fecun*H*(1-(L/Ka))-(dev_L+del_L)*L
    dA <-  dev_L*L - (omeg*Hu + del_A)*A
    dH <-  omeg*Hu*A - (gon + del_A)*H
    # return the rate of change
    list(c(dL, dA, dH))
  }) # end with(as.list ...
}

likelihood <- function(x) # forzamientos para el solver de la ode
{ 
  if(x[1] <= 0 ){
    print("Negative param")
    res = -86829146000
  }else{
    # print("Positive param")
    pars <- c(fecun = fec,Ka = K,Hu = Hum,omeg = x[1],
              del_L = delta_L,del_A = delta_A,dev_L = d_L,
              gon = a) # death rate group 2
    
    sd_t <- x[2]
    
    population <- c(L=10.0, A=0.0,H= 0.0) #Vector inicial para ODE
    
    # z <- ode(y=population,
    # times = 0:nrow(y), func = "derivs", method = "ode45",
    # parms = parms, dllname = "model_vec_test1",
    # initfunc = "initmod", nout = 1,
    # outnames = "Sum")
    z <- ode(y = population, times = 0:nrow(y), func = vect, parms = pars)
    
    #Aquí corre el ODE
    
    colnames(z)[2:4] <- c("L", "A", "Ah")
    
    z <- as.data.frame(z)
    z <- z[-1, ]
    
    L <- y$L
    A <- y$A
    Ah <- y$Ah
    
    res <-  sum(dnorm(L, mean = z$L, sd = sd_t, log = T))  +
      sum(dnorm(A, mean = z$A , sd = sd_t, log = T))  +
      sum(dnorm(Ah, mean = z$Ah, sd = sd_t, log = T))  #cálculo de la loglikelihood en función de las desviaciones estándar
    # sum(dnorm((A+Ah), mean = (z$A + z$Ah), sd = sd, log = T))  
  }
  
  # print("res:")
  # print(res)
  return(res)
}

#------------------------------------------------------------------------------#
# Params values:
fec = 200
K = 250000
Hum = 1600000
omega_t = 0.2
delta_L = 0.2
delta_A = 0.25
d_L = 0.08
a = 1/8

feasability_cond(delta_L,delta_A,d_L,a,fec ,K,Hum,omega_t)
vec_eq <- eq_point(delta_L, delta_A, d_L, a, fec, K, Hum, omega_t)
stab_cond = function(x){return(feasability_cond(delta_L,delta_A,d_L,a,fec ,K,Hum,x))}
x <- seq(0.1, 1, by=.05)
eq_mat <- matrix(0, ncol = 3, nrow = length(x))
for(i in c(1:length(x))){
  eq_mat[i,] <- eq_point(delta_L, delta_A, d_L, a, fec, K, Hum, x[i])
}
# Stability in function of omega.
x <- seq(0, 100, by=.05)
vec_stab <- lapply(x, stab_cond)
stab_cond_vec <- data.frame(omega = x, stab = as.array(vec_stab, nrow =1, ncol = length(x)))
plot(stab_cond_vec)

# Equilibrium points as a function of omega.
eq_point_df <- data.frame(omega = x, L_eq = eq_mat[,1], A_eq = eq_mat[,2], Ah_eq = eq_mat[,3])
df_plot2 <- reshape2::melt(eq_point_df, id.vars = c("omega"))
df_plot_fil <- df_plot2 %>% filter(variable == "L_eq" | variable == "Ah_eq" )
df_plot_fil <- df_plot2 %>% filter(variable == "A_eq"  )
ggplot(df_plot_fil,aes(omega, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Number of individuals") +
  ggtitle("Equilibrium points")+
  scale_color_manual(name = "",
                     labels = c("Adult mosquito",  "Adult handling mosquito"),
                     values=c('#FF00F6','#2F822B'))+
  theme_bw()

x <- seq(0, 100, by=.05)
omeg_coeff <- lapply(x, coeff_omega)
coeff_omega <- data.frame(omega = x, stab = as.array(omeg_coeff, nrow =1, ncol = length(x)))
plot(coeff_omega)


######------------------RHO Data Analysis---------------------######
Path = "~/MAD_MODEL/VECTOR_MODEL/data/df_rho.dat"
df_rho <- data.frame(t(read.table(Path, header=FALSE)))
colnames(df_rho) <- c("time", "rho", "date")
df_rho$date = as.Date(df_rho$date , "%Y-%m-%d")
init_date <- df_rho$date[1]
end_date <-  tail(df_rho$date, n=1)
n_days = as.numeric(end_date - init_date) 
df_date <- data.frame(date = seq(init_date, end_date, "days"), time = seq(0,n_days,1))

df_rho$time <- NULL
df_rho$rho <- as.numeric(df_rho$rho)
ggplot(df_rho) + 
  geom_line(aes(x = date, y =rho)) +
  ggtitle("Encounter rate computed from Citizen Science model") +
  theme_bw()
avg_rho <- mean(df_rho$rho[lubridate::year(df_rho$date) == 2018 ])

######-------------BG COUNTS ESTIMATE ---------------#####
#Upload data from bgtraps.
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
  ylab("Adults mosquitoes estimated in BCN") +
  theme_bw()

# Moving Average:
l = length(bg_traps$prop)
dt = 7
mov = mov_avg(l,dt,bg_traps$value)
df_bg_mov <- data.frame(date = bg_traps$date, mov_bg_count = mov)

l = length(df_rho$rho)
dt = 7
mov = mov_avg(l,dt,df_rho$rho)

df_rho_mov <- data.frame(date = df_rho$date, rho = mov)
# df_rho_mov <- merge(df_rho_mov, df_date, by = "time")
df_tot <- merge(df_rho_mov, df_bg_mov, by = "date")

# Scatter plot with marginal distributions.
# Total BG counts:
ggscatterstats(data = df_tot,
               x = mov_bg_count,
               y = rho, 
               xlab ="Adult mosquito per m^2 (BG traps)" , 
               ylab =expression(rho),
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

linear_mod <- lm(df_tot$rho~df_tot$mov_bg_count)
summary(linear_mod)

# Total counts/population size (past approach):
df_tot$density <- df_tot$mov_bg_count/pop_density
linear_mod <- lm(df_tot$rho~df_tot$density)
summary(linear_mod)

ggscatterstats(data = df_tot,
               x = density,
               y = rho, 
               xlab ="Adult mosquito per m^2 (BG traps)" , 
               ylab =expression(rho),
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 


# Set to zero (artificially) the number of mosquitoes in winter:
head(df_tot)
df_tot$mov_bg_count[lubridate::month(df_tot$date) == 1 ] <- 0
df_tot$mov_bg_count[lubridate::month(df_tot$date) == 2 ] <- 0
df_tot$mov_bg_count[lubridate::month(df_tot$date) == 3 ] <- 0
df_tot$mov_bg_count[lubridate::month(df_tot$date) == 12 ] <- 0
ggscatterstats(data = df_tot,
               x = mov_bg_count,
               y = rho, 
               xlab ="Adult mosquito per m^2 (BG traps)" , 
               ylab =expression(rho),
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

linear_mod <- lm(df_tot$rho~df_tot$mov_bg_count)
summary(linear_mod)

# Check different trends in different years:
head(df_tot)
df_tot_fil <- df_tot %>% filter(lubridate::year(date) == 2018) 

ggscatterstats(data = df_tot_fil,
               x = mov_bg_count,
               y = rho, 
               xlab ="Adult mosquito per m^2 (BG traps)" , 
               ylab =expression(rho),
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

linear_mod <- lm(df_tot_fil$rho~df_tot_fil$mov_bg_count)
summary(linear_mod)

head(df_tot)
df_tot_fil <- df_tot %>% filter(lubridate::year(date) == 2019) 

ggscatterstats(data = df_tot_fil,
               x = mov_bg_count,
               y = rho, 
               xlab ="Adult mosquito per m^2 (BG traps)" , 
               ylab =expression(rho),
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

linear_mod <- lm(df_tot_fil$rho~df_tot_fil$mov_bg_count)
summary(linear_mod)

head(df_tot)
df_tot_fil <- df_tot %>% filter(lubridate::year(date) == 2020) 

ggscatterstats(data = df_tot_fil,
               x = mov_bg_count,
               y = rho, 
               xlab ="Adult mosquito per m^2 (BG traps)" , 
               ylab =expression(rho),
               ggplot.component = list(ggplot2::geom_rug(sides = "b"))) 

linear_mod <- lm(df_tot_fil$rho~df_tot_fil$mov_bg_count)
summary(linear_mod)

# Multiple regression one plot:
df_tot$year <- as.character(lubridate::year(df_tot$date))
ggplot(df_tot, aes(mov_bg_count, rho, shape=year, colour=year, fill=year)) +
  geom_smooth(method="lm") +
  geom_point(size=3) +
  theme_bw() + 
  xlab("Estimation of adult mosquitoes in BCN") +
  ylab("rho(A)")

#---------------INTEGRATION----------------#
times <- seq(0, 100, by = 0.01)
parameters <- c(fecun = fec,Ka = K,Hu = Hum,omeg = omega_t,del_L = delta_L,del_A = delta_A,dev_L = d_L,gon = a)
# state <- c(L = vec_eq[1], A = vec_eq[2], H = vec_eq[3])
state <- c(L = 10, A = -0, H = 0)
out2 <- as.data.frame(ode(y = state, times = times, func = vect, parms = parameters))
df_plot2 <- reshape2::melt(out2, id.vars = c("time"))
ggplot(df_plot2,aes(time, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Counts") +
  ggtitle("Vector dynamics")+
  scale_color_manual(name = "",
                     labels = c("Larva", "Adult mosquito", "Adult handling mosquito"),
                     values=c('#FF00F6','#FF2C00','#2F822B'))+
  theme_bw()


# Pseudo Data to check the oprimization method.
slopevalues = function(x){return(likelihood(c(x, trueSD)))}
slopevalues(0.1)
x <- seq(-1, 10, by=.05)
slopelikelihoods = lapply(x, slopevalues )
plot (x, slopelikelihoods , type="l", xlab = "values of slope parameter a", ylab = "Log likelihood")