rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")

# Create Pseudo Data:
gam1 = 0.2
sd = 1
x <- seq(1,100,0.1)
l <- length(x)
fun <- function(gam1){
  return(gam1*x + 3)
}
y <- gam1*x + 3 + rnorm(l,mean =0, sd = sd)
# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE

true1 = 0.2
# Likelihood function
likelihood <- function(par) # forzamientos para el solver de la ode
{ 
  z <- fun(par[1])
  sd <- par[2]
  res <- #cálculo de la loglikelihood en función de las desviaciones estándar
    sum(dnorm(y, mean = z, sd = sd, log = T)) 
  
  # print(paste0("res:", res))
  return(res)
}

# Example: plot the likelihood profile of the slope.
slopevalues <- function(par){
  return(likelihood(c(par,sd)))
} 

vec <- seq(0, 1, by=.01)
slopelikelihoods <- lapply(vec, slopevalues )

plot(vec, slopelikelihoods , type="l", xlab = "values of slope gamma 1", ylab = "Log likelihood")

# Prior distribution
prior = function(param){
  a = param[1]
  sd = param[2]
  aprior = dunif(a, min=0, max=4, log = T)
  sdprior = dunif(sd, min=0, max=10, log = T)
  return(aprior+sdprior)
}

# Posterior distribution (sum because we work with logarithms)
posterior = function(param){
  return (likelihood(y) + prior(param))
}


######## Metropolis algorithm ################

proposalfunction = function(param){
  vec <- c(rnorm(1, mean = param, sd= 0.1)
           ,abs(rnorm(1,mean = param[],sd = 0.3)))
  return(vec)
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,2))
  chain[1,] = startvalue
  prop_mat <- vector("numeric", length = iterations)
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    # print("Proposal:")
    # print(proposal)
    # print(paste0("Iteration:",i))
    if(is.na(posterior(proposal))){
      print("na likelyhood")      
    }
    probab = exp(posterior(proposal) - posterior(chain[i,]))
    # print("posterior(proposal,ob_data,forcs_mat):",posterior(proposal,ob_data,forcs_mat))
    prop_mat[i] <- probab
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,] } 
  }
  return(chain)
}

startvalue = c(0.1,0.21)
iterations = 10000
chain = run_metropolis_MCMC(startvalue, iterations)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))


### Summary: #######################

par(mfrow = c(2,4))
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]))
abline(v = true1, col="red" )

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = true1, col="red" )¡

# for comparison:
summary(lm(y~x))
