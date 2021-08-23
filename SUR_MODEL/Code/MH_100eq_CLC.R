rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")

# Create Pseudo data:
# Create Pseudo Data:
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model.c")
dyn.load("model.so")

gam1 = 0.2
gam2 = 1.2
gam3 = 3
# We create a vector with the constant parameters.
parms = c(gam1,gam2,gam3)
# We set the initial conditions to zero.
Y <- c(y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0, y5 = 0, y6 = 0, y7 = 0, y8 = 0, y9 = 0, y10 = 0,
       y11 = 0.0, y12 = 0.0, y13 = 0.0, y14 = 0, y15 = 0, y16 = 0, y17 = 0, y18 = 0, y19 = 0, y20 = 0,
       y21 = 0.0, y22 = 0.0, y23 = 0.0, y24 = 0, y25 = 0, y26 = 0, y27 = 0, y28 = 0, y29 = 0, y30 = 0,
       y31 = 0.0, y32 = 0.0, y33 = 0.0, y34 = 0, y35 = 0, y36 = 0, y37 = 0, y38 = 0, y39 = 0, y40 = 0,
       y41 = 0.0, y42 = 0.0, y43 = 0.0, y44 = 0, y45 = 0, y46 = 0, y47 = 0, y48 = 0, y49 = 0, y50 = 0,
       y51 = 0.0, y52 = 0.0, y53 = 0.0, y54 = 0, y55 = 0, y56 = 0, y57 = 0, y58 = 0, y59 = 0, y60 = 0,
       y61 = 0.0, y62 = 0.0, y63 = 0.0, y64 = 0, y65 = 0, y66 = 0, y67 = 0, y68 = 0, y69 = 0, y70 = 0,
       y71 = 0.0, y72 = 0.0, y73 = 0.0, y74 = 0, y75 = 0, y76 = 0, y77 = 0, y78 = 0, y79 = 0, y80 = 0,
       y81 = 0.0, y82 = 0.0, y83 = 0.0, y84 = 0, y85 = 0, y86 = 0, y87 = 0, y88 = 0, y89 = 0, y90 = 0,
       y91 = 0.0, y92 = 0.0, y93 = 0.0, y94 = 0, y95 = 0, y96 = 0, y97 = 0, y98 = 0, y99 = 0, y100 = 0)

Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")

# List with the data frames of the forcings, sort as the c code.
forcs_mat <- list(data.matrix(down))

min_t <- min(down$time)
max_t <- max(down$time)
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model",
           initfunc = "initmod", nout = 1,
           outnames = "Sum", initforc = "forcc",
           forcings = down, 
           fcontrol = list(method = "constant")) 


ode <- data.frame(out)
ode$Sum <- NULL
head(ode)

saveRDS(ode, file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo_5eq.rds")
# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model.c")
dyn.load("model.so")

true1 = 0.2
true2 = 1.2
true3 = 3
trueSD = 1
# Likelihood function
likelihood <- function(y, #datos
                       x, # vector con los parámetros
                       forcings) # forzamientos para el solver de la ode
{ 
  
  pars <- c(gam1 = x[1], # death rate group 1
            gam2 = x[2],
            gam3 = x[3]) # death rate group 2
  
  sd <- x[4]
  
  population <- c(y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0.0, y5 = 0.0) #Vector inicial para ODE
  
  forcs_mat <- list(data.matrix(forcings))
  
  z <- ode(y = population,
           times = 0:nrow(y), func = "derivs", method = "ode45",
           dllname = "model", initfunc = "initmod", nout = 0, 
           parms = pars, initforc = "forcc", forcings = forcs_mat, 
           fcontrol = list(method = "constant")) #Aquí corre el ODE
  
  colnames(z)[2:6] <- c("P1", "P2", "P3", "P4", "P5")
  
  z <- as.data.frame(z)
  z <- z[-1, ]
  
  P1 <- y$X1
  P2 <- y$X2
  P3 <- y$X3
  P4 <- y$X4
  P5 <- y$X5
  
  res <- #cálculo de la loglikelihood en función de las desviaciones estándar
    sum(dnorm(P1, mean = z$P1, sd = sd, log = T)) +
    sum(dnorm(P2, mean = z$P2, sd = sd, log = T)) +
    sum(dnorm(P3, mean = z$P3, sd = sd, log = T))+
    sum(dnorm(P4, mean = z$P4, sd = sd, log = T)) +
    sum(dnorm(P5, mean = z$P5, sd = sd, log = T))
  
  return(res)
}


# Registrations:
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")
head(down)
forcs_mat <- data.matrix(down)

# Pseudo Data to check the oprimization method.
ob_data <- readRDS(file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo_5eq.rds")
colnames(ob_data) <- c("time", "X1", "X2", "X3", "X4", "X5")
l <- nrow(ob_data)
ob_data$X1 <- ob_data$X1 + rnorm(l,0,trueSD)
ob_data$X2 <- ob_data$X2 + rnorm(l,0,trueSD)
ob_data$X3 <- ob_data$X3 + rnorm(l,0,trueSD)
ob_data$X4 <- ob_data$X4 + rnorm(l,0,trueSD)
ob_data$X5 <- ob_data$X5 + rnorm(l,0,trueSD)

head(ob_data)
# Example: plot the likelihood profile of the slope.
slopevalues <- function(par){
  return(likelihood(ob_data,c(par, true2,true3, trueSD),forcs_mat))
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
  aprior = dnorm(a, sd=1,  log = T)
  bprior = dnorm(b, sd=1,  log = T)
  cprior = dnorm(c, sd=1,  log = T)
  sdprior = dnorm(sd, sd=1,  log = T)
  return(aprior+bprior+cprior+sdprior)
}

# Posterior distribution (sum because we work with logarithms)
posterior = function(param, y, forc){
  return (likelihood(y,param,forc) + prior(param))
}


######## Metropolis algorithm ################

proposalfunction = function(param){
  return(abs(rnorm(4, mean = param, sd= c(0.1,0.5,0.3,0.4))))
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,4))
  chain[1,] = startvalue
  prop_mat <- vector("numeric", length = iterations)
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    #print("Proposal:")
    #print(proposal)
    probab = exp(posterior(proposal,ob_data,forcs_mat) - posterior(chain[i,],ob_data,forcs_mat))
    prop_mat[i] <- probab
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,] } 
  }
  return(chain)
}

startvalue = c(0.1,0.1,0.1,0.1)
iterations = 10000
chain = run_metropolis_MCMC(startvalue, iterations)

burnIn = 5000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
### Summary: #######################

par(mfrow = c(2,4))
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]))
abline(v = true1, col="red" )
hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]))
abline(v = true2, col="red" )
hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),3]) )
abline(v = trueSD, col="red" )

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = true1, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = true2, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSD, col="red" )

filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/chain_MH_",iterations,".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain, file = filename)

print("Optimization finish")

