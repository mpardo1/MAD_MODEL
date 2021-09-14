rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("coda")
# Create Pseudo data:
# Create Pseudo Data:
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_4eq.c")
dyn.load("model_4eq.so")

gam1 = 0.2
gam2 = 1.2
gam3 = 3
# We create a vector with the constant parameters.
parms = c(gam1,gam2,gam3)
# We set the initial conditions to zero.
Y <- c(y1 = 0, y2=0, y3=0, y4 = 0)
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")

# List with the data frames of the forcings, sort as the c code.
forcs_mat <- list(data.matrix(down))

min_t <- min(down$time)
max_t <- max(down$time)
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model_4eq",
           initfunc = "initmod", nout = 1,
           outnames = "Sum", initforc = "forcc",
           forcings = down, 
           fcontrol = list(method = "constant")) 


ode <- data.frame(out)
ode$Sum <- NULL
head(ode)

saveRDS(ode, file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo_4eq.rds")
# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE
true1 = gam1
true2 = gam2
true3 = gam3
trueSD = 1
# Likelihood function
likelihood <- function(x) # forzamientos para el solver de la ode
{ 
  # print("x:")
  # print(x)
  if(x[1] < 0 | x[2] < 0 | x[3] < 0 ){
    print("Negative param")
    res = -86829146000
  }else{
    pars <- c(gam1 = x[1], # death rate group 1
              gam2 = x[2],
              gam3 = x[3]) # death rate group 2
    
    sd <- x[4]
    
    population <- c(y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0.0) #Vector inicial para ODE
    
    forcs_mat <- list(data.matrix(forcings))
    
    z <- ode(y = population,
             times = 0:nrow(y), func = "derivs", method = "ode45",
             dllname = "model_4eq", initfunc = "initmod", nout = 0, 
             parms = pars, initforc = "forcc", forcings = forcs_mat, 
             fcontrol = list(method = "constant")) #Aquí corre el ODE
    
    colnames(z)[2:5] <- c("P1", "P2", "P3", "P4")
    
    z <- as.data.frame(z)
    z <- z[-1, ]
    
    P1 <- y$X1
    P2 <- y$X2
    P3 <- y$X3
    P4 <- y$X4
    
    res <- #cálculo de la loglikelihood en función de las desviaciones estándar
      sum(dnorm(P1, mean = z$P1, sd = sd, log = T)) +
      sum(dnorm(P2, mean = z$P2, sd = sd, log = T)) +
      sum(dnorm(P3, mean = z$P3, sd = sd, log = T)) +
      sum(dnorm(P4, mean = z$P4, sd = sd, log = T)) 
  }
  
  # print("res:")
  # print(res)
  return(res)
}


# Registrations:
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")
head(down)
forcs_mat <- data.matrix(down)

# Pseudo Data to check the oprimization method.
ob_data <- readRDS(file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo_4eq.rds")
colnames(ob_data) <- c("time", "X1", "X2", "X3", "X4")
l <- nrow(ob_data)
ob_data$X1 <- ob_data$X1 + rnorm(l,0,trueSD)
ob_data$X2 <- ob_data$X2 + rnorm(l,0,trueSD)
ob_data$X3 <- ob_data$X3 + rnorm(l,0,trueSD)
ob_data$X4 <- ob_data$X4 + rnorm(l,0,trueSD)

y <- ob_data
forcings <- down
head(ob_data)

# Prior distribution
prior = function(param){
  a = param[1]
  b = param[2]
  c = param[3]
  sd = param[4]
  aprior = dnorm(a, sd=1,  log = T)
  bprior = dnorm(b, sd=1,  log = T)
  cprior = dnorm(c, sd=1,  log = T)
  sdprior = dnorm(sd, sd=1,  log = T)
  return(aprior+bprior+cprior+sdprior)
}

######## Metropolis algorithm ################

proposalfunction = function(param){
  vec <- param + c(rnorm(3, mean = c(0,0,0), sd= c(0.1,0.1,0.1))
                   ,abs(rnorm(1,mean = 0 ,sd = 0.3)))
  return(vec)
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,4))
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


startvalue = c(0.1,1,2.5,0.5)
iterations = 100000
chain = run_metropolis_MCMC(startvalue, iterations)

filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/chain_MH_op_4eq_3param",iterations,".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain, file = filename)

# chain <- mcmc(chain)
# summary(chain)
# plot(chain)
# 
# library(BayesianTools)
# correlationPlot(data.frame(chain))
# 
# # Convergence diagnosis:
# 
# print("Optimization finish")
chain2 = run_metropolis_MCMC(startvalue, 10000)
filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/chain2_MH_op_4eq_3param",iterations,".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain2, file = filename)
# combinedchains = mcmc.list(chain, chain2)
# plot(combinedchains)
# gelman.diag(combinedchains)
# gelman.plot(combinedchains)
