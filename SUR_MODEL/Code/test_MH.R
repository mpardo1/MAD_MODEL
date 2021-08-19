rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")

# Create Pseudo Data:
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_test.c")
dyn.load("model_test.so")

gam1 = 0.2
# We create a vector with the constant parameters.
parms = c(gam1)
# We set the initial conditions to cero.
Y <- c(y1 = 0)
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")

# List with the data frames of the forcings, sort as the c code.
forcs_mat <- list(data.matrix(down))

min_t <- min(down$time)
max_t <- max(down$time)
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model_test",
           initfunc = "initmod", nout = 1,
           outnames = "Sum", initforc = "forcc",
           forcings = down, 
           fcontrol = list(method = "constant")) 


ode <- data.frame(out)
ode$Sum <- NULL
head(ode)

saveRDS(ode, file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo.rds")
# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_test.c")
dyn.load("model_test.so")

true1 = 0.2
true2 = 0.01
true3 = 1.4
trueSD = 1
# Likelihood function
likelihood <- function(y, #datos
                       x, # vector con los parámetros
                       forcings) # forzamientos para el solver de la ode
{ 
  
  pars <- c(gam1 = x[1]) # death rate group 3
  
  sd <- x[2]
  
  population <- c(y1 = 0.0) #Vector inicial para ODE
  
  forcs_mat <- list(data.matrix(forcings))
  
  z <- ode(y = population,
           times = 0:nrow(y), func = "derivs", method = "ode45",
           dllname = "model_test", initfunc = "initmod", nout = 0, 
           parms = pars, initforc = "forcc", forcings = forcs_mat, 
           fcontrol = list(method = "constant")) #Aquí corre el ODE
  
  colnames(z)[2] <- c("P1")
  
  z <- as.data.frame(z)
  z <- z[-1, ]
  
  P1 <- y$X1
  res <- #cálculo de la loglikelihood en función de las desviaciones estándar
    sum(dnorm(P1, mean = z$P1, sd = sd, log = T)) 
  
  # print(paste0("res:", res))
  return(res)
}



# Registrations:
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")
head(down)
forcs_mat <- data.matrix(down)

# Pseudo Data to check the oprimization method.
ob_data <- readRDS(file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo.rds")
ob_data <- ob_data[1:300,1:2]
colnames(ob_data) <- c("time", "X1")
l <- nrow(ob_data)
ob_data$X1 <- ob_data$X1 + rnorm(l,0,trueSD)


head(ob_data)
# Example: plot the likelihood profile of the slope.
slopevalues <- function(par){
  return(likelihood(ob_data,c(par,trueSD),forcs_mat))
} 

vec <- seq(0, 1, by=.01)
slopelikelihoods <- lapply(vec, slopevalues )

plot(vec, slopelikelihoods , type="l", xlab = "values of slope gamma 1", ylab = "Log likelihood")

# Prior distribution
prior = function(param){
  a = param[1]
  sd = param[2]
  aprior = dnorm(a, sd=5,  log = T)
  sdprior = dunif(sd, min=0, max=10, log = T)
  return(aprior+sdprior)
}

# Posterior distribution (sum because we work with logarithms)
posterior = function(param, y, forc){
  return (likelihood(y,param,forc) + prior(param))
}


######## Metropolis algorithm ################

proposalfunction = function(param){
  vec <- param + rnorm(2, mean = 0, sd = 0.1)
  return(vec)
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  prop_mat <- vector("numeric", length = iterations)
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    print("Proposal:")
    print(proposal)
    # print("Proposal:")
    # print(proposal)
    # print(paste0("Iteration:",i))
    if(is.na(posterior(proposal,ob_data,forcs_mat))){
      print("na likelyhood")      
    }
    probab = exp(posterior(proposal,ob_data,forcs_mat) - posterior(chain[i,],ob_data,forcs_mat))
    # print("posterior(proposal,ob_data,forcs_mat):",posterior(proposal,ob_data,forcs_mat))
    
    print("probab:")
    print(probab)
    prop_mat[i] <- probab
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,] 
      } 
  }
  return(chain)
}

startvalue = c(0,0.21)
iterations = 10000
chain = run_metropolis_MCMC(startvalue, iterations)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))


### Summary: #######################

# par(mfrow = c(2,4))
# hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of a", xlab="True value = red line" )
# abline(v = mean(chain[-(1:burnIn),1]))
# abline(v = true1, col="red" )
# hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
# abline(v = mean(chain[-(1:burnIn),2]))
# abline(v = true, col="red" )
# 
# plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
# abline(h = true1, col="red" )
# plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
# abline(h = trueSD, col="red" )

# for comparison:
# summary(lm(y~x))

filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/chain_MH_",iterations,".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain, file = filename)

print("Optimization finish")