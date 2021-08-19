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


# Función que calcula la loglikelihood del modelo ----------------------------------

# Esta será la función a optimizar

ll_ode <- function(x, # vector con los parámetros
                   forcings, # forzamientos para el solver de la ode
                   y, # datos
                   devs){ #desviaciones estándar para calcular la loglikelihood
  
  pars <- c(gam1 = x[1]) # death rate group 3
  
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
    sum(dnorm(z$P1 - P1, sd = devs[1], log = T))
}

# Carga datos -------------------------------------------------------------

# Registrations:
Path = "~/MAD_MODEL/SUR_MODEL/data/Observed_data_2300.data"
# ob_data <- data.frame(t(read.table(Path, header=FALSE)))

ob_data <- readRDS(file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo.rds")
colnames(ob_data) <- c("time", "X1")
head(ob_data)

# f_inicio <- as.Date("2021-05-14")
# f_fin <- as.Date("2021-06-12")
# 
n_inicio <- min(ob_data$time)
n_fin <- max(ob_data$time)

input2 <- ob_data[n_inicio:n_fin, ]


# Estima desviación estandar de cada serie --------------------------------

# Esto nos permite asignar una likelihood al conjunto del modelo, considerando
# las series independientes.

devs <- c()
spl <- (input2$X1)
fit <- smooth.spline(x = 1:nrow(input2), y = spl, df = 4)
devs[1] <- sd(spl - predict(fit)$y)

# Forzamientos ------------------------------------------------------------

# Registrations:
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")
head(down)

# Hospitalizados iniciales ------------------------------------------------

# A punto de empezar el proceso -------------------------------------------

best <- -999999999 #LL inicial a mejorar

# load("seeds_CAN.RData") #Cargamos seeds de los valores iniciales de los pars.
seeds <- matrix( runif(100,0,1), ncol = 100, nrow = 1)
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

Cores <- 19 #Numero de cores a utilizar.
it <- 0
while(condition){
  #Ahora viene la paralelización
  parall <- mclapply(1:sims, mc.cores = Cores, mc.preschedule = F,function(k){
    it <- it + 1
    
    fit <- optim(par = seeds[, k], fn = ll_ode, forcings = down, y = input2, 
                 devs = devs,lower=c(0), upper=rep(Inf, 1), control = list(fnscale = -1, maxit = 500, parscale = seeds[, k]))
    
    if((k %% 1000) == 0) {
      cat("This was seed no. ", k, "\n")
      cat("This fit: ", fit$value, "\n")
    }
    
    fit
  })
  
  lhs <- parall
  
  rm(parall) #Para evitar fugas de memoria
  
  filename <- paste0("param_MAD_MODEL", round, ".RData") #Salva cada ronda de optimizaciones, por si acaso
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
  parmat <- matrix(NA, nrow = length(index), ncol = 1)
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
