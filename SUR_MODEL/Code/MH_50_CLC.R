rm(list = ls())
start_time <- Sys.time()

library("parallel")
library("tidyverse")
library("deSolve")

# Create Pseudo data:
# Create Pseudo Data:
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_50eq.c")
dyn.load("model_50eq.so")

gam1 = 0.2
gam2 = 1.2
gam3 = 3
# We create a vector with the constant parameters.
parms = c(gam1,gam2,gam3)
# We set the initial conditions to zero.
Y <-  c(y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0, y5 = 0, y6 = 0, y7 = 0, y8 = 0, y9 = 0, y10 = 0,
        y11 = 0.0, y12 = 0.0, y13 = 0.0, y14 = 0, y15 = 0, y16 = 0, y17 = 0, y18 = 0, y19 = 0, y20 = 0,
        y21 = 0.0, y22 = 0.0, y23 = 0.0, y24 = 0, y25 = 0, y26 = 0, y27 = 0, y28 = 0, y29 = 0, y30 = 0,
        y31 = 0.0, y32 = 0.0, y33 = 0.0, y34 = 0, y35 = 0, y36 = 0, y37 = 0, y38 = 0, y39 = 0, y40 = 0,
        y41 = 0.0, y42 = 0.0, y43 = 0.0, y44 = 0, y45 = 0, y46 = 0, y47 = 0, y48 = 0, y49 = 0, y50 = 0)
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")

# List with the data frames of the forcings, sort as the c code.
forcs_mat <- list(data.matrix(down))

min_t <- min(down$time)
max_t <- max(down$time)
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model_50eq",
           initfunc = "initmod", nout = 1,
           outnames = "Sum", initforc = "forcc",
           forcings = down, 
           fcontrol = list(method = "constant")) 


ode <- data.frame(out)
ode$Sum <- NULL
head(ode)

saveRDS(ode, file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo_50eq.rds")
# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_50eq.c")
dyn.load("model_50eq.so")

true1 = 0.2
true2 = 1.2
true3 = 3
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
    
    population <- c(y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0, y5 = 0, y6 = 0, y7 = 0, y8 = 0, y9 = 0, y10 = 0,
                    y11 = 0.0, y12 = 0.0, y13 = 0.0, y14 = 0, y15 = 0, y16 = 0, y17 = 0, y18 = 0, y19 = 0, y20 = 0,
                    y21 = 0.0, y22 = 0.0, y23 = 0.0, y24 = 0, y25 = 0, y26 = 0, y27 = 0, y28 = 0, y29 = 0, y30 = 0,
                    y31 = 0.0, y32 = 0.0, y33 = 0.0, y34 = 0, y35 = 0, y36 = 0, y37 = 0, y38 = 0, y39 = 0, y40 = 0,
                    y41 = 0.0, y42 = 0.0, y43 = 0.0, y44 = 0, y45 = 0, y46 = 0, y47 = 0, y48 = 0, y49 = 0, y50 = 0)
    
    forcs_mat <- list(data.matrix(forcings))
    
    z <- ode(y = population,
             times = 0:nrow(y), func = "derivs", method = "ode45",
             dllname = "model_50eq", initfunc = "initmod", nout = 0, 
             parms = pars, initforc = "forcc", forcings = forcs_mat, 
             fcontrol = list(method = "constant")) #Aquí corre el ODE
    
    colnames(z)[2:51] <- c("P1", "P2", "P3", "P4", "P5","P6", "P7", "P8", "P9", "P10",
                           "P11", "P12", "P13", "P14", "P15","P16", "P17", "P18", "P19", "P20",
                           "P21", "P22", "P23", "P24", "P25","P26", "P27", "P28", "P29", "P30",
                           "P31", "P32", "P33", "P34", "P35","P36", "P37", "P38", "P39", "P40",
                           "P41", "P42", "P43", "P44", "P45","P46", "P47", "P48", "P49", "P50")
    
    z <- as.data.frame(z)
    z <- z[-1, ]
    
    P1 <- y$X1;  P2 <- y$X2;  P3 <- y$X3;  P4 <-  y$X4; P5 <- y$X5;  P6 <- y$X6  ;P7 <- y$X7;  P8 <- y$X8;   P9 <- y$X9;P10 <- y$X10
    P11 <- y$X11;P12 <- y$X12;P13 <- y$X13;P14 <- y$X14;P15 <- y$X15;P16 <- y$X16;P17 <- y$X17;P18 <- y$X18; P19 <- y$X19;P20 <- y$X20
    P21 <- y$X21;P22 <- y$X22;P23 <- y$X23;P24 <- y$X24;P25 <- y$X25;P26 <- y$X26;P27 <- y$X27;P28 <- y$X28; P29 <- y$X29;P30 <- y$X30
    P31 <- y$X31;P32 <- y$X32;P33 <- y$X33;P34 <- y$X34;P35 <- y$X35;P36 <- y$X36;P37 <- y$X37;P38 <- y$X38; P39 <- y$X39;P40 <- y$X40
    P41 <- y$X41;P42 <- y$X42;P43 <- y$X43;P44 <- y$X44;P45 <- y$X45;P46 <- y$X46;P47 <- y$X47;P48 <- y$X48; P49 <- y$X49;P50 <- y$X50
    
    res <- #cálculo de la loglikelihood en función de las desviaciones estándar
      sum(dnorm(P1, mean = z$P1, sd = sd, log = T)) + sum(dnorm(P6, mean = z$P6, sd = sd, log = T)) +
      sum(dnorm(P2, mean = z$P2, sd = sd, log = T)) + sum(dnorm(P7, mean = z$P7, sd = sd, log = T)) +
      sum(dnorm(P3, mean = z$P3, sd = sd, log = T)) + sum(dnorm(P8, mean = z$P8, sd = sd, log = T)) +
      sum(dnorm(P4, mean = z$P4, sd = sd, log = T)) + sum(dnorm(P9, mean = z$P9, sd = sd, log = T)) +
      sum(dnorm(P5, mean = z$P5, sd = sd, log = T)) + sum(dnorm(P10, mean = z$P10, sd = sd, log = T)) +
      sum(dnorm(P11, mean = z$P11, sd = sd, log = T)) + sum(dnorm(P16, mean = z$P16, sd = sd, log = T)) +
      sum(dnorm(P12, mean = z$P12, sd = sd, log = T)) + sum(dnorm(P17, mean = z$P17, sd = sd, log = T)) +
      sum(dnorm(P13, mean = z$P13, sd = sd, log = T)) + sum(dnorm(P18, mean = z$P18, sd = sd, log = T)) +
      sum(dnorm(P14, mean = z$P14, sd = sd, log = T)) + sum(dnorm(P19, mean = z$P19, sd = sd, log = T)) +
      sum(dnorm(P15, mean = z$P15, sd = sd, log = T)) + sum(dnorm(P20, mean = z$P20, sd = sd, log = T)) +
      sum(dnorm(P21, mean = z$P21, sd = sd, log = T)) + sum(dnorm(P26, mean = z$P26, sd = sd, log = T)) +
      sum(dnorm(P22, mean = z$P22, sd = sd, log = T)) + sum(dnorm(P27, mean = z$P27, sd = sd, log = T)) +
      sum(dnorm(P23, mean = z$P23, sd = sd, log = T)) + sum(dnorm(P28, mean = z$P28, sd = sd, log = T)) +
      sum(dnorm(P24, mean = z$P24, sd = sd, log = T)) + sum(dnorm(P29, mean = z$P29, sd = sd, log = T)) +
      sum(dnorm(P25, mean = z$P25, sd = sd, log = T)) + sum(dnorm(P30, mean = z$P30, sd = sd, log = T)) +
      sum(dnorm(P31, mean = z$P31, sd = sd, log = T)) + sum(dnorm(P36, mean = z$P36, sd = sd, log = T)) +
      sum(dnorm(P32, mean = z$P32, sd = sd, log = T)) + sum(dnorm(P37, mean = z$P37, sd = sd, log = T)) +
      sum(dnorm(P33, mean = z$P33, sd = sd, log = T)) + sum(dnorm(P38, mean = z$P38, sd = sd, log = T)) +
      sum(dnorm(P34, mean = z$P34, sd = sd, log = T)) + sum(dnorm(P39, mean = z$P39, sd = sd, log = T)) +
      sum(dnorm(P35, mean = z$P35, sd = sd, log = T)) + sum(dnorm(P40, mean = z$P40, sd = sd, log = T)) +
      sum(dnorm(P41, mean = z$P41, sd = sd, log = T)) + sum(dnorm(P46, mean = z$P46, sd = sd, log = T)) +
      sum(dnorm(P42, mean = z$P42, sd = sd, log = T)) + sum(dnorm(P47, mean = z$P47, sd = sd, log = T)) +
      sum(dnorm(P43, mean = z$P43, sd = sd, log = T)) + sum(dnorm(P48, mean = z$P48, sd = sd, log = T)) +
      sum(dnorm(P44, mean = z$P44, sd = sd, log = T)) + sum(dnorm(P49, mean = z$P49, sd = sd, log = T)) +
      sum(dnorm(P45, mean = z$P45, sd = sd, log = T)) + sum(dnorm(P50, mean = z$P50, sd = sd, log = T)) 
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
ob_data <- readRDS(file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo_50eq.rds")
colnames(ob_data) <- c("time", "X1", "X2", "X3", "X4", "X5" , "X6", "X7", "X8", "X9", "X10",
                       "X11", "X12", "X13", "X14", "X15" , "X16", "X17", "X18", "X19", "X20",
                       "X21", "X22", "X23", "X24", "X25" , "X26", "X27", "X28", "X29", "X30",
                       "X31", "X32", "X33", "X34", "X35" , "X36", "X37", "X38", "X39", "X40",
                       "X41", "X42", "X43", "X44", "X45" , "X46", "X47", "X48", "X49", "X50")
l <- nrow(ob_data)
ob_data[,2:51] <- ob_data[,2:51] + matrix(rnorm(l*50,0,trueSD), ncol = 50, nrow = l)


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
  vec <- param + c(rnorm(3, mean = c(0,0,0), sd= c(0.00001,0.00001,0.00005))
                   ,abs(rnorm(1,mean = 0 ,sd = 0.00003)))
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
    # print("Iteration:")
    # print(i)
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

startvalue = c(0.15,1.2,2.7,0.5)
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

filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/chain_MH_50eq_3param_",iterations,"_",Sys.Date(),".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain, file = filename)

print("Optimization finish")
end_time <- Sys.time()
diff_time <- end_time - start_time
print("Execution time:")
print(diff_time)
