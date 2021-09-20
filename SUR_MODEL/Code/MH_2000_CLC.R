rm(list = ls())
start_time <- Sys.time()

library("parallel")
library("tidyverse")
library("deSolve")

# Create Pseudo data:
# Create Pseudo Data:
# Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# # Path = paste(PC,Path, sep="")
# 
# setwd(Path)
# system("R CMD SHLIB model_2000eq.c")
# dyn.load("model_2000eq.so")
# 
# gam1 = 0.2
# gam2 = 1.2
# gam3 = 3
# # We create a vector with the constant parameters.
# parms = c(gam1,gam2,gam3)
# # We set the initial conditions to zero.
# Y <-  matrix(0, nrow = 1, ncol=2000)
# Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
# down <- data.frame(t(read.table(Path, header=FALSE)))
# colnames(down) <- c("time", "down")
# 
# # List with the data frames of the forcings, sort as the c code.
# forcs_mat <- list(data.matrix(down))
# 
# min_t <- min(down$time)
# max_t <- max(down$time)
# times <- seq(min_t,max_t, 1)
# 
# out <- ode(Y, times, func = "derivs",
#            parms = parms, dllname = "model_2000eq",
#            initfunc = "initmod", nout = 1,
#            outnames = "Sum", initforc = "forcc",
#            forcings = down, 
#            fcontrol = list(method = "constant")) 
# 
# 
# ode <- data.frame(out)
# ode$Sum <- NULL
# head(ode)
# ode_mat <- as.matrix(ode)
# saveRDS(ode, file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo_2000eq.rds")

# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_2000eq.c")
dyn.load("model_2000eq.so")

true1 = 0.0250265
true2 = 0.146925
true3 = 0.482396
trueSD = 1
# Likelihood function
likelihood <- function(y, #datos
                       x, # vector con los parámetros
                       forcings) # forzamientos para el solver de la ode
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
    
    population <- matrix(0,nrow = 1,ncol=2000)
    
    forcs_mat <- list(data.matrix(forcings))
    
    z <- ode(y = population,
             times = 0:nrow(y), func = "derivs", method = "ode45",
             dllname = "model_2000eq", initfunc = "initmod", nout = 0, 
             parms = pars, initforc = "forcc", forcings = forcs_mat, 
             fcontrol = list(method = "constant")) #Aquí corre el ODE
    
    
    z <- as.data.frame(z)
    z <- z[-1, ]
    
    res <- 0
    for(i in c(1:2000)){
      res <- res +  sum(dnorm( y[,i+1], mean = z[,i+1], sd = sd, log = T))
    }
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

# Read Observed data:
ob_data <-read.table("~/MAD_MODEL/SUR_MODEL/Code/Observed_data_2300.data", header=FALSE, sep= " ")
ob_data <- t(ob_data[,1:2000])
# 
# # Pseudo Data to check the oprimization method.
# ob_data <- readRDS(file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo_2000eq.rds")
# l <- nrow(ob_data)
# # ob_data[,2:1001] <- ob_data[,2:1001] + matrix(rnorm(l*2000,0,trueSD), ncol = 2000, nrow = l)
# 
# 
# head(ob_data)
# # Example: plot the likelihood profile of the slope.
# slopevalues <- function(par){
#   return(likelihood(ob_data,c(par, true2,true3, trueSD),forcs_mat))
# } 
# 
# vec <- seq(0, 1, by=.1)
# slopelikelihoods <- lapply(vec, slopevalues )
# 
# plot(vec, slopelikelihoods , type="l", xlab = "values of slope gamma 1", ylab = "Log likelihood")

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
  vec <- param + c(rnorm(3, mean = c(0,0,0), sd= c(0.00008,0.00008,0.00008))
                   ,abs(rnorm(1,mean = 0 ,sd = 0.00008)))
  return(vec)
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,4))
  like = array(dim = c(iterations+1,1))
  chain[1,] = startvalue
  prop_mat <- vector("numeric", length = iterations)
  like[1] <- posterior(chain[1,],ob_data,forcs_mat)
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    print("Iteration:")
    print(i)
    like1 <- posterior(proposal,ob_data,forcs_mat)
    like2 <- like[i]
    probab = exp(like1 - like2)
    prop_mat[i] <- probab
    if (runif(1) < probab){
      chain[i+1,] = proposal
      like[i+1] = like1
    }else{
      chain[i+1,] = chain[i,] 
      like[i+1] = like2
    } 
  }
  return(chain)
}

startvalue = c(0.02,0.1,0.4,0.5)
iterations = 10000
chain = run_metropolis_MCMC(startvalue, iterations)

burnIn = 50
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
### Summary: #######################
# 
# par(mfrow = c(2,4))
# hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of a", xlab="True value = red line" )
# abline(v = mean(chain[-(1:burnIn),1]))
# abline(v = true1, col="red" )
# hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of b", xlab="True value = red line")
# abline(v = mean(chain[-(1:burnIn),2]))
# abline(v = true2, col="red" )
# hist(chain[-(1:burnIn),3],nclass=30, main="Posterior of c", xlab="True value = red line")
# abline(v = mean(chain[-(1:burnIn),3]))
# abline(v = true3, col="red" )
# hist(chain[-(1:burnIn),4],nclass=30, main="Posterior of sd", xlab="True value = red line")
# abline(v = mean(chain[-(1:burnIn),4]) )
# abline(v = trueSD, col="red" )
# 
# plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
# abline(h = true1, col="red" )
# plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
# abline(h = true2, col="red" )
# plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of c", )
# abline(h = true3, col="red" )
# plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
# abline(h = trueSD, col="red" )
# 
filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/chain_MH_2000eq_3param_",iterations,"_",Sys.Date(),".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain, file = filename)

print("Optimization finish")
end_time <- Sys.time()
diff_time <- end_time - start_time
print("Execution time:")
print(diff_time)
