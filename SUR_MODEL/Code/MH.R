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

true1 = 0.02
true2 = 0.02
true3 = 0.03
trueSD = 10
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
  P2 <- y$X2
  P3 <- y$X3
  
  res <- #cálculo de la loglikelihood en función de las desviaciones estándar
    sum(dnorm(P1, mean = z$P1, sd = sd, log = T)) +
    sum(dnorm(P2, mean = z$P2, sd = sd, log = T)) +
    sum(dnorm(P3, mean = z$P3, sd = sd, log = T))
}


# Registrations:
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")
head(down)
forcs_mat <- data.matrix(down)

# Example: plot the likelihood profile of the slope.
slopelikelihoods = likelihood(mat,par,forcs_mat)
plot (seq(3, 7, by=.05), slopelikelihoods , type="l", xlab = "values of slope parameter a", ylab = "Log likelihood")


# Prior distribution
prior = function(param){
  a = param[1]
  b = param[2]
  c = param[3]
  sd = param[4]
  aprior = dunif(a, min=0, max=4, log = T)
  bprior = dunif(b, min=0, max=4, log = T)
  cprior = dunif(c, min=0, max=4, log = T)
  sdprior = dunif(sd, min=0, max=30, log = T)
  return(aprior+bprior+cprior+sdprior)
}
