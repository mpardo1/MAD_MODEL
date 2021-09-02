rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("coda")

###############   ODE INTEGRATION   ##################
Path = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_vec_cte.c")
dyn.load("model_vec_cte.so")

f = 200
K = 250000
H = 160000
omega_t = 4
trueSD = 1
# We create a vector with the constant parameters.
parms = c(f,K,H, omega_t)
# We set the initial conditions to cero.
Y <- c(y1 = 100.0, y2 = 0.0, y3 = 0.0)
min_t <- 1
max_t <- 365
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model_vec_cte",
           initfunc = "initmod", nout = 1,
           outnames = "Sum") 

ode <- data.frame(out) 
ode$Sum <- NULL
saveRDS(ode, file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo.rds")


likelihood <- function(x) # forzamientos para el solver de la ode
{ 
  if(x[1] <= 0 ){
    print("Negative param")
    res = -86829146000
  }else{
    pars <- c(f = f,K = K,H = H,omega = x[1]) # death rate group 2
    
    sd <- x[2]
    
    population <- c(y1 = 100.0, y2 = 0.0, y3 = 0.0) #Vector inicial para ODE
    
    
    z <- ode(y=population,
             times = 0:nrow(y), func = "derivs", 
             dllname = "model_vec_cte" , parms = pars,
             initfunc = "initmod") 
    #Aquí corre el ODE
    
    colnames(z)[2:4] <- c("L", "A", "Ah")
    
    z <- as.data.frame(z)
    z <- z[-1, ]
    
    L <- y$L
    A <- y$A
    Ah <- y$Ah
    
    res <-  sum(dnorm(L, mean = z$L, sd = 1, log = T))  +
            sum(dnorm(A, mean = z$A , sd = 1, log = T))  +
            sum(dnorm(Ah, mean = z$Ah, sd = 1, log = T))  #cálculo de la loglikelihood en función de las desviaciones estándar
      # sum(dnorm((A+Ah), mean = (z$A + z$Ah), sd = sd, log = T))  
  }
  
  # print("res:")
  # print(res)
  return(res)
}

# Pseudo Data to check the oprimization method.
ob_data <- readRDS(file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo.rds")
colnames(ob_data) <- c("time", "L", "A", "Ah")
l <- nrow(ob_data)
ob_data$L <- ob_data$L + rnorm(l,0,trueSD)
ob_data$A <- ob_data$A + rnorm(l,0,trueSD)
ob_data$Ah <- ob_data$Ah + rnorm(l,0,trueSD)

true1 = omega_t

y <- ob_data

slopevalues = function(x){return(likelihood(c(x, trueSD)))}
x <- seq(0.5, 100, by=.05)
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

