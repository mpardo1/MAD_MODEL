rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("coda")

# Params values:
fec = 20
K = 43
Hum = 13
omega_t = 0.2
delta_L = 0.5
delta_A = 3
d_L = 0.8
a = 0.01

# Equilibrium points:
eq_point <- function(delta_L, delta_A, d_L, a, fec, K, Hum, omega_t){
  Ah_eq <-  K*(((omega_t*Hum*d_L)/((a+delta_A)*(omega_t*Hum +delta_A)))-((1/(a*fec))*(d_L+delta_L)))
  L_eq <- ((a*fec*Ah_eq)/((a*fec/K)*Ah_eq+(d_L+delta_L)))
  A_eq <- ((d_L/(omega_t*Hum+delta_A))*((a*fec*Ah_eq)/(((a*fec/K)*Ah_eq)+(d_L+delta_L))))
  eq <- c(L_eq, A_eq, Ah_eq)
}


vec_eq <- eq_point(delta_L, delta_A, d_L, a, fec, K, Hum, omega_t)

ode_value <- function(L, A, Ah){
  val <- c()
  val[1] = a*fec*Ah*(1-(L/K))-(d_L+delta_L)*L
  val[2] = d_L*L - (omega_t*Hum + delta_A)*A
  val[3] = omega_t*Hum*A - (a + delta_A)*Ah
  return(val)
}

ode_value(vec_eq[1],vec_eq[2],vec_eq[3])
ode_value(0,0,0)

# Feasability condition:
feasability_cond <- function(del_L, del_A, dev_L, gon, fec, Kar, Hum, omega){
  fes <- FALSE
  val <- (gon*fec*omega*Hum*dev_L)/((gon + del_A)*(omega*Hum + del_A)*(dev_L + del_L))
  if(val > 1){
    fes <- TRUE
  }
  return(fes)
}

feasability_cond(delta_L,delta_A,d_L,a,fec ,K,Hum,omega_t)
###############   ODE INTEGRATION   ##################
Path = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_vec_test.c")
dyn.load("model_vec_test.so")


trueSD = 1
# We create a vector with the constant parameters.
parms = c(fec = fec,Ka = K,H = Hum,omeg = omega_t,del_L = delta_L,del_A = delta_A,dev_L = d_L,gon = a)
# We set the initial conditions to zero.
Y <- c(y1 = L_eq, y2 = A_eq, y3 = Ah_eq)
min_t <- 1
max_t <- 4200
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs", method = "ode45",
           parms = parms, dllname = "model_vec_test",
           initfunc = "initmod", nout = 1,
           outnames = "Sum") 

ode <- data.frame(out) 
ode$Sum <- NULL


# saveRDS(ode, file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo_cte.rds")
df_plot <- reshape2::melt(ode, id.vars = c("time"))
ggplot(df_plot,aes(time, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Counts") +
  ggtitle("Vector dynamics")+
  scale_color_manual(name = "",
                     labels = c("Larva", "Adult mosquito", "Adult handling mosquito"),
                     values=c('#FF00F6','#FF2C00','#2F822B'))+
  theme_bw()

#-----------------------------------------------------------------------------------
# ODE system in R:
vect <- function(t, state, parameters) {
   with(as.list(c(state, parameters)),{
     # rate of change
     
       dL <-  gon*fecun *H*(1-(L/Ka))-(dev_L+del_L)*L
       dA <-  dev_L*L - (omeg*H + del_A)*A
       dH <-  omeg*H*A - (gon + del_A)*H
         # return the rate of change
         list(c(dL, dA, dH))
       }) # end with(as.list ...
}

fec = 20
K = 43
Hum = 13
omega_t = 0.2
delta_L = 0.5
delta_A = 3
d_L = 0.8
a = 0.01

vec_eq <- eq_point(delta_L, delta_A, d_L, a, fec, K, Hum, omega_t)

times <- seq(0, 100, by = 0.01)
parameters <- c(fecun = fec,Ka = K,H = Hum,omeg = omega_t,del_L = delta_L,del_A = delta_A,dev_L = d_L,gon = a)
state <- c(L = vec_eq[1], A = vec_eq[2], H = vec_eq[3])
# state <- c(L = -2222.025, A = -317.4321, H = 13)
out <- as.data.frame(ode(y = state, times = times, func = vect, parms = parameters))
df_plot <- reshape2::melt(out, id.vars = c("time"))
ggplot(df_plot,aes(time, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Counts") +
  ggtitle("Vector dynamics")+
  scale_color_manual(name = "",
                     labels = c("Larva", "Adult mosquito", "Adult handling mosquito"),
                     values=c('#FF00F6','#FF2C00','#2F822B'))+
  theme_bw()

#-----------------------------------------------------------------------------------
likelihood <- function(x) # forzamientos para el solver de la ode
{ 
  if(x[1] <= 0 ){
    print("Negative param")
    res = -86829146000
  }else{
    pars <- c(f = f,K = K,H = H,omega = x[1]) # death rate group 2
    
    sd_t <- x[2]
    
    population <- c(y1 = 10.0, y2 = 0.0, y3 = 0.0) #Vector inicial para ODE
    
    z <- ode(y=population,
             times = 0:nrow(y), func = "derivs", method = "ode45",
             parms = parms, dllname = "model_vec_test",
             initfunc = "initmod", nout = 1,
             outnames = "Sum") 

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

# Pseudo Data to check the oprimization method.
ob_data <- readRDS(file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo_cte.rds")
colnames(ob_data) <- c("time", "L", "A", "Ah")
l <- nrow(ob_data)
# ob_data$L <- ob_data$L + rnorm(l,0,trueSD)
# ob_data$A <- ob_data$A + rnorm(l,0,trueSD)
# ob_data$Ah <- ob_data$Ah + rnorm(l,0,trueSD)

true1 = omega_t

y <- ob_data

slopevalues = function(x){return(likelihood(c(x, trueSD)))}
slopevalues(0.1)
x <- seq(0.5, 10, by=.05)
slopelikelihoods = lapply(x, slopevalues )
plot (x, slopelikelihoods , type="l", xlab = "values of slope parameter a", ylab = "Log likelihood")
# Prior distribution
prior = function(param){
  a = param[1]
  sd = param[2]
  aprior = dnorm(a, sd=1,  log = T)
  sdprior = dnorm(sd, sd=1,  log = T)
  return(aprior+sdprior)
}

######## Metropolis algorithm ################

proposalfunction = function(param){
  vec <- c(rnorm(1, mean = param[1], sd= 0.1)
           ,abs(rnorm(1,mean = param[2] ,sd = 0.1)))
  
  return(vec)
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    print("Iteration:")
    print(i)
    probab = exp(likelihood(proposal)+ prior(proposal) - likelihood(chain[i,])- prior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,] } 
  }
  return(chain)
}


startvalue = c(0.03,0.7)
iterations = 10000
chain = run_metropolis_MCMC(startvalue, iterations)

filename <- paste0("~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/chain_MH_op",iterations,".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain, file = filename)
# burnIn = 500
# acceptance1 = 1-mean(duplicated(chain[-(1:burnIn),]))

chain <- mcmc(chain)
# summary(chain)
# plot(chain)

# library(BayesianTools)
# correlationPlot(data.frame(chain))
#
# # Convergence diagnosis:
#
# print("Optimization finish")
chain2 = run_metropolis_MCMC(startvalue, iterations)
# burnIn = 5000
# acceptance2 = 1-mean(duplicated(chain2[-(1:burnIn),]))
# 
# chain2 <- mcmc(chain2)
filename <- paste0("~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/chain2_MH_op",iterations,".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain, file = filename)

# 
# combinedchains = mcmc.list(chain, chain2)
# plot(combinedchains)
# gelman.diag(combinedchains)
# gelman.plot(combinedchains)

