rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("coda")

# Params values:
fec = 2000
K = 43000
Hum = 13554
omega_t = 0.2
delta_L = 0.2
delta_A = 0.3
d_L = 8
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
system("R CMD SHLIB model_vec_test1.c")
dyn.load("model_vec_test1.so")

trueSD = 1
# We create a vector with the constant parameters.
parms = c(fecun = fec, Ka = K, Hu = Hum, omeg = omega_t, del_L = delta_L, del_A = delta_A, dev_L = d_L, gon = a)
# We set the initial conditions to zero.
eps = 0.001
veq_eq <- vec_eq + eps
Y <- c(y1 = vec_eq[1], y2 = vec_eq[2], y3 = vec_eq[3])
# Y <- c(y1 = 10, y2 = 0, y3 = 0)
min_t <- 1
max_t <- 100
times <- seq(min_t,max_t, 1)
out1 <- as.data.frame(ode(Y, times, func = "derivs", method = "ode45",
           parms = parms, dllname = "model_vec_test1",
           initfunc = "initmod", nout = 1,
           outnames = "Sum"))

# out1$y1 <- NULL
# out1$y2 <- NULL
out1$Sum <- NULL

df_plot <- reshape2::melt(out1, id.vars = c("time"))
ggplot(df_plot,aes(time, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Counts") +
  ggtitle("Vector dynamics")+
  scale_color_manual(name = "",
                     labels = c("Larva", "Adult mosquito", "Adult handling mosquito"),
                     values=c('#FF00F6','#FF2C00','#2F822B')) +
  theme_bw()

#-----------------------------------------------------------------------------------
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


# Params values:
fec = 2000
K = 4300
Hum = 13554
omega_t = 0.2
delta_L = 0.2
delta_A = 0.3
d_L = 14
a = 0.01
feasability_cond(delta_L,delta_A,d_L,a,fec ,K,Hum,omega_t)

vec_eq <- eq_point(delta_L, delta_A, d_L, a, fec, K, Hum, omega_t)

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
saveRDS(out2, file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo_cte.rds")
#-----------------------------------------------------------------------------------
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
x <- seq(-1, 10, by=.05)
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

