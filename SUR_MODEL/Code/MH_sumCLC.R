rm(list = ls())
start_time <- Sys.time()

library("parallel")
library("tidyverse")
library("deSolve")

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
  if(x[1] < 0 | x[2] < 0 | x[3] < 0 | x[4] == 0){
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
    
    # Sum the columns of the ode integration:
    sum_group1 <- rowSums(z[,2:31])
    sum_group2 <- rowSums(z[,31:648])
    sum_group3 <- rowSums(z[,648:2000])
    
    # Remove the df of solutions z:
    rm(z)
    
    # Compute the LL:
      res <- sum(dnorm( y$sum1, mean = sum_group1, sd = sd, log = T)) + 
              sum(dnorm( y$sum2, mean = sum_group2, sd = sd, log = T)) + 
              sum(dnorm( y$sum3, mean = sum_group3, sd = sd, log = T)) 
    
  }
  
  print("res:")
  print(res)
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
sum_1 <- rowSums(ob_data[,2:31])
sum_2 <- rowSums(ob_data[,31:648])
sum_3 <- rowSums(ob_data[,648:2000])

ob_data <- data.frame(time = ob_data[,1], sum1 = sum_1 , sum2 = sum_2 , sum3 = sum_3 )
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
  vec <- param + c(rnorm(3, mean = c(0,0,0), sd= c(0.08,0.08,0.08))
                   ,abs(rnorm(1,mean = 0 ,sd = 0.0008)))
  return(vec)
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,4))
  like = array(dim = c(iterations+1,1))
  chain[1,] = startvalue
  prop_mat <- vector("numeric", length = iterations)
  like[1] <- posterior(chain[1,],ob_data,forcs_mat)
  # print("chain:")
  # print(chain[1,])
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    print("Iteration:")
    print(i)
    like1 <- posterior(proposal,ob_data,forcs_mat)
    like2 <- like[i]
    probab = exp(like1 - like2)
    prop_mat[i] <- probab
    # print("probab:")
    # print(probab)
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

startvalue = c(0,0,0,1)
iterations = 15000
chain = run_metropolis_MCMC(startvalue, iterations)

burnIn = 50

acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
print(paste0("Acceptance rate: ", acceptance))

filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/chain_sum_2000eq_3param_",iterations,"_",Sys.Date(),".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain, file = filename)

print("Optimization finish")
end_time <- Sys.time()
diff_time <- end_time - start_time
print("Execution time:")
print(diff_time)
