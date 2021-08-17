rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model.c")
dyn.load("model.so")

true1 = 0.2
true2 = 0.01
true3 = 1.4
trueSD = 1
# Likelihood function
likelihood <- function(y, #datos
                       x, # vector con los parámetros
                       forcings) # forzamientos para el solver de la ode
                       { 
  
  pars <- c(gam1 = x[1], # death rate group 1
            gam2 = x[2], # death rate group 2
            gam3 = x[3]) # death rate group 3
  
  sd <- x[4]
  
  population <- c(y1 = 0.0, y2 = 0.0, y3 = 0.0) #Vector inicial para ODE
  
  forcs_mat <- list(data.matrix(forcings))
 
  z <- ode(y = population,
           times = 0:nrow(y), func = "derivs", method = "ode45",
           dllname = "model", initfunc = "initmod", nout = 0, 
           parms = pars, initforc = "forcc", forcings = forcs_mat, 
           fcontrol = list(method = "constant")) #Aquí corre el ODE
  
  colnames(z)[2:4] <- c("P1", "P2", "P3")
  
  z <- as.data.frame(z)
  z <- z[-1, ]
  
  P1 <- y$X1
  print(paste0("P1", P1))
  P2 <- y$X2
  P3 <- y$X3
  
  res <- #cálculo de la loglikelihood en función de las desviaciones estándar
    sum(dnorm(P1, mean = z$P1, sd = sd, log = T)) +
    sum(dnorm(P2, mean = z$P2, sd = sd, log = T)) +
    sum(dnorm(P3, mean = z$P3, sd = sd, log = T))
 
  return(res)
}


# Registrations:
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")
head(down)
forcs_mat <- data.matrix(down)

# Pseudo Data to check the oprimization method.
ob_data <- readRDS(file = "ode_pseudo.rds")
colnames(ob_data) <- c("time", "X1", "X2", "X3")
l <- nrow(ob_data)
ob_data$X1 <- ob_data$X1 + rnorm(l,0,trueSD)
ob_data$X2 <- ob_data$X2 + rnorm(l,0,trueSD)
ob_data$X3 <- ob_data$X3 + rnorm(l,0,trueSD)


head(ob_data)
# Example: plot the likelihood profile of the slope.
slopevalues <- function(par){
  return(likelihood(ob_data,c(par, true2,true3,trueSD),forcs_mat))
} 

vec <- seq(0, 1, by=.01)
slopelikelihoods <- lapply(vec, slopevalues )

plot(vec, slopelikelihoods , type="l", xlab = "values of slope gamma 1", ylab = "Log likelihood")

# Prior distribution
prior = function(param){
  a = param[1]
  b = param[2]
  c = param[3]
  sd = param[4]
  aprior = dunif(a, min=0, max=4, log = T)
  bprior = dunif(b, min=0, max=4, log = T)
  cprior = dunif(c, min=0, max=4, log = T)
  sdprior = dunif(sd, min=0, max=10, log = T)
  return(aprior+bprior+cprior+sdprior)
}

# Posterior distribution (sum because we work with logarithms)
posterior = function(param, y, forc){
  return (likelihood(y,param,forc) + prior(param))
}


######## Metropolis algorithm ################

proposalfunction = function(param){
  return(abs(rnorm(4, mean = param, sd= c(0.1,0.5,0.3))))
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,4))
  chain[1,] = startvalue
  prop_mat <- vector("numeric", length = iterations)
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    print("Proposal:")
    print(proposal)
    probab = exp(posterior(proposal,ob_data,forcs_mat) - posterior(chain[i,],ob_data,forcs_mat))
    prop_mat[i] <- probab
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,] } 
    }
  return(list(chain,prop_mat))
}

startvalue = c(0.1,0.21,1,0.1)
iterations = 100000
chain = run_metropolis_MCMC(startvalue, iterations)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))

filename <- "~/MAD_MODEL/SUR_MODEL/Code/chain_MH.RData" #Salva cada ronda de optimizaciones, por si acaso
save(chain, file = filename)


