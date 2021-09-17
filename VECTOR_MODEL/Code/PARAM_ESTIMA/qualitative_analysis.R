rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("coda")

#-------------------------------FUNCTIONS--------------------------------------#
# Equilibrium points:
eq_point <- function(delta_L, delta_A, d_L, a, fec, K, Hum, omega_t){
  Ah_eq <-  K*(((omega_t*Hum*d_L)/((a+delta_A)*(omega_t*Hum +delta_A)))-((1/(a*fec))*(d_L+delta_L)))
  L_eq <- ((a*fec*Ah_eq)/((a*fec/K)*Ah_eq+(d_L+delta_L)))
  A_eq <- ((d_L/(omega_t*Hum+delta_A))*((a*fec*Ah_eq)/(((a*fec/K)*Ah_eq)+(d_L+delta_L))))
  eq <- c(L_eq, A_eq, Ah_eq)
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
x <- seq(0, 100, by=.05)
eq_mat <- as.matrix(0, ncol = 3, nrow = length(x))
for(i in range(1:length(x))){
  eq_mat[i] <- eq_point(delta_L, delta_A, d_L, a, fec, K, Hum, x[i])
}
# Stability in function of omega.
x <- seq(0, 100, by=.05)
vec_stab <- lapply(x, stab_cond)
stab_cond_vec <- data.frame(x <- x, stab = as.array(vec_stab, nrow =1, ncol = length(x)))
plot(stab_cond_vec)

# Equilibrium points as a function of omega.
eq_point_vec <- lapply(x, eq_point_ap)
stab <- as.array(eq_point_vec)
eq_point_df <- data.frame(omega = x, L_eq = stab[])

plot(eq_point_df)
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