rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")

# Load C code:
Path = "~/MAD_MODEL/SUR_MODEL/Code/PARAM_ESTIMATION/MH/"
setwd(Path)
system("R CMD SHLIB model_2000eq.c")
dyn.load("model_2000eq.so")

# Real value parameters:
gam1 = 0.0250265
gam2 = 0.146925
gam3 = 0.482396

# Función que calcula la loglikelihood del modelo ----------------------------------
# Esta será la función a optimizar

ll_ode <- function(x, # vector con los parámetros
                   forcings, # forzamientos para el solver de la ode
                   y, # datos
                   devs){ #desviaciones estándar para calcular la loglikelihood
  
  if(x[1] < 0 | x[2] < 0 | x[3] < 0){
    res = -86829146000
  }else{
    pars <- c(gam1 = x[1], # death rate group 1
              gam2 = x[2], # death rate group 2
              gam3 = x[3]) # death rate group 3
    
    population <- matrix(0,nrow = 1,ncol=2000) #Vector inicial para ODE
    
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
    res <- sum(dnorm( y[,1], mean = sum_group1, sd = devs[1], log = T)) + 
      sum(dnorm( y[,2], mean = sum_group2, sd = devs[2], log = T)) + 
      sum(dnorm( y[,3], mean = sum_group3, sd = devs[3], log = T)) 
    
  }
  return(res)
}

# Carga datos -------------------------------------------------------------
# Registrations:
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")
head(down)
forcs_mat <- data.matrix(down)

# Read Observed data:
ob_data <-read.table("~/MAD_MODEL/SUR_MODEL/data/Observed_data_2300.data", header=FALSE, sep= " ")
ob_data <- t(ob_data[,1:2000])

# n_inicio <- which(data$FECHA == f_inicio)
# n_fin <- which(data$FECHA == f_fin)

n_inicio <- 1
n_fin <- 2000

input2 <- matrix(0, nrow = nrow(ob_data), ncol = 3)


# Estima desviación estandar de cada serie --------------------------------

# Esto nos permite asignar una likelihood al conjunto del modelo, considerando
# las series independientes.
input2[,1] <- rowSums(ob_data[,2:31])
input2[,2] <- rowSums(ob_data[,31:648])
input2[,3] <- rowSums(ob_data[,648:2000])

devs <- c()
spl <- input2[,1]
fit <- smooth.spline(x = 1:length(spl), y = spl, df = 4)
devs[1] <- sd(spl - predict(fit)$y)

spl <- input2[,2]
fit <- smooth.spline(x = 1:length(spl), y = spl, df = 4)
devs[2] <- sd(spl - predict(fit)$y)

spl <- input2[,3]
fit <- smooth.spline(x = 1:length(spl), y = spl, df = 4)
devs[3] <- sd(spl - predict(fit)$y)
for(i in c(1:3)){
  if(devs[i] == 0){
    devs[i] <- 0.001
  }
}

# A punto de empezar el proceso -------------------------------------------

best <- -999999999 #LL inicial a mejorar

# load("seeds_CAN.RData") #Cargamos seeds de los valores iniciales de los pars.
seeds <- matrix( runif(900,0,1), ncol = 300, nrow = 3)
sols <- NA #Pre-aloco el número de combinaciones paramétricas en 2 unidades de LL de la mejor

set.seed(476468713)

condition <- T 
#Esta condición se hará falsa en el bucle de después, si en las rondas la LL no
#mejora más de 2 puntos y hay más de 1000 combinaciones paramétricas a menos de
#2 puntos de loglikelihood de la mejor solución.

round <- 1
# Paralelización

# sims <- 2 #Número de combinaciones paramétricas a explorar
sims <- ncol(seeds) #Número de combinaciones paramétricas a explorar

Cores <- 1
# Cores <- parallel::detectCores()#Numero de cores a utilizar.
it <- 0
while(condition){
  start_time <- Sys.time()
  print(paste0("Iteration: ", it))
  #Ahora viene la paralelización
  parall <- mclapply(1:sims, mc.cores = Cores, mc.preschedule = F,function(k){
    
    
    fit <- optim(par = seeds[, k], fn = ll_ode, forcings = down, y = input2, 
                 devs = devs, control = list(fnscale = -1, maxit = 500, parscale = seeds[, k]))
    
    if((k %% 1000) == 0) {
      cat("This was seed no. ", k, "\n")
      cat("This fit: ", fit$value, "\n")
    }
    
    fit
  })
  end_time <- Sys.time()
  diff_time <- end_time - start_time
  print("Execution time iteration i:")
  print(diff_time)
  it <- it + 1
  lhs <- parall
  
  rm(parall) #Para evitar fugas de memoria
  
  filename <- paste0("~/MAD_MODEL/SUR_MODEL/OUTPUT/NM/param_MAD_MODEL_1core_900it", round, ".RData") #Salva cada ronda de optimizaciones, por si acaso
  save(lhs, file = filename)
  
  # Ahora, recuperamos la loglikelihood de cada combinación de parámetros
  logl <- rep(NA, sims)
  for(i in 1:sims) logl[i] <- lhs[[i]]$value
  
  # Evaluamos las condiciones para parar el bucle
  best2 <- max(logl, na.rm = T)
  pc1 <- best < best2
  pc2 <- best > (best2 - 2)
  
  sols <- sum(logl > (max(logl, na.rm = T) - 2), na.rm = T)
  pc3 <- sols > 1000
  
  condition <- pc1 * pc2 * pc3
  condition <- !condition
  
  # Seleccionamos las mejores combinaciones de parámetros para mandar una nueva
  # ronda, cogemos las combinaciones que estén a 2 unidades de distancia de la
  # mejor, o en su defecto, las 250 mejores combinaciones.
  if(sols < 250){
    index <- order(logl, decreasing = T)[1:250]
  } else {
    index <- order(logl, decreasing = T)[1:sols]
  }
  
  n <- 1
  parmat <- matrix(NA, nrow = length(index), ncol = 3)
  for(i in index){
    parmat[n, ] <- lhs[[i]]$par
    n <- n + 1
  }
  
  # Muestreamos las mejores combinaciones de parámetros, cada uno de ellos
  # independientemente. De esta forma estamos rompiendo las posibles
  # correlaciones entre los parámetros, y en algún sentido, hacemos la
  # aproximación menos bayesiana al no samplear la distribución multidimensional
  # de los parámetros.
  seeds <- t(apply(X = parmat, MARGIN = 2, sample, size = 100000, replace = T))
  
  rm(parmat)
  rm(lhs)
  
  print(paste0("Round: ", round, ";  Best likelihood: ", best2, 
               "; Solutions: ", sols, "\n"))
  
  best <- best2
  round <- round + 1
} 
print("The optimization finished.")
