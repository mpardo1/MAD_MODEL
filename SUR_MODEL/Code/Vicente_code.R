
library(deSolve)
library(parallel)
library(tidyverse)

# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE
dyn.load("pcrind.so")


# Función que calcula la loglikelihood del modelo ----------------------------------

# Esta será la función a optimizar

ll_ode <- function(x, # vector con los parámetros
                   forcings, # forzamientos para el solver de la ode
                   y, # datos
                   devs, #desviaciones estándar para calcular la loglikelihood
                   hosp){ #hospitalizados totales y críticos
  
  pars <- c(b = x[1], #infectividad per capita
            fi = x[2], #infectividad relativa de los presintomáticos
            eI = x[3], #contención de los poco sintomáticos
            eY = x[4], #contención de los hospitalizados
            sigma = x[5], #tasa de incubación
            gamma1 = x[6], #tasa aparición síntomas
            p = x[7], # probabilidad de hospitalización
            kappa = x[8], #detectabilidad
            gamma2 = x[9], #tasa recuperación asintomáticos
            alpha = x[10], #tasa de entrada en uci
            delta = x[11], #tasa de mortalidad
            gamma3 = x[16], #tasa de recuperación hospitalizados
            gamma4 = x[17], #tasa de recuperación críticos
            deltab = x[18]) #tasa de mortalidad basal
  
  # Población Cantabria: 582905
  susc <- 582905 - sum(x[12:15]) - hosp[1] - x[19]
  
  population <- c(susc, x[12], x[13], x[14], x[15], hosp[1] - hosp[2], 
                  hosp[2], 0, x[19], 0, 0, 0) #Vector inicial para ODE
  
  z <- ode(y = population,
           times = 0:nrow(y), func = "derivs", method = "ode45",
           dllname = "pcrind", initfunc = "initmod", nout = 0, 
           parms = pars, initforc = "forcc", forcings = forcings, 
           fcontrol = list(method = "constant")) #Aquí corre el ODE
  
  colnames(z)[2:13] <- c("S", "E", "I1", "A", "Ad", "I2", "Y", "D_ac", "R",
                         "Y_ac", "I2_ac", "Pos_ac")
  
  z <- as.data.frame(z)
  z <- z[-1, ]
  
  D_ac <- (y$exitus)
  I2 <- y$ingressats
  Y_ac <- (y$novas_uci)
  Y <- y$critics
  I2_ac <- (y$ingressos)
  Pos_ac <- (y$positius)
  
  res <- #cálculo de la loglikelihood en función de las desviaciones estándar
    sum(dnorm(diff(c(0, z$D_ac)) - D_ac, sd = devs[1], log = T)) +
    sum(dnorm((z$I2 + z$Y) - I2, sd = devs[2], log = T)) +
    sum(dnorm(diff(c(0, z$Y_ac)) - Y_ac, sd = devs[3], log = T)) +
    sum(dnorm(z$Y - Y, sd = devs[4], log = T)) +
    sum(dnorm(diff(c(0, z$I2_ac)) - I2_ac, sd = devs[5], log = T)) +
    sum(dnorm(diff(c(0, z$Pos_ac)) - Pos_ac, sd = devs[6], log = T))
  
  #Penalizaciones parámetros
  penaltyb <- ((x > 1e-5) * 1.0  - 1) * 10^6
  
  penalty2 <- ((x[2] > 1) * -1.0) * 10^6
  penalty3 <- ((x[3] > 1) * -1.0) * 10^6
  penalty4 <- ((x[4] > 1) * -1.0) * 10^6
  penalty7 <- ((x[7] > .6) * -1.0) * 10^6
  penalty12 <- ((x[12] < 10) * -1.0) * 10^6
  penalty13 <- ((x[13] < 10) * -1.0) * 10^6
  penalty14 <- ((x[14] < 10) * -1.0) * 10^6
  penalty15 <- ((x[15] < 10) * -1.0) * 10^6
  p12 <- ((x[12] > (347.2431 * 1.2)) * -1.0) * 10^6
  p13 <- ((x[13] > (167.59675 * 1.2)) * -1.0) * 10^6
  p14 <- ((x[14] > (792.1791 * 1.2)) * -1.0) * 10^6
  p15 <- ((x[15] > (726.7072 * 1.2)) * -1.0) * 10^6
  
  
  p2 <- ((x[5] > 2) * -1.0) * 10^6
  p3 <- ((x[6] > 2) * -1.0) * 10^6
  p4 <- ((x[8] > 1100) * -1.0) * 10^6
  p5 <- ((x[9] > 2) * -1.0) * 10^6
  p6 <- ((x[10] > 2) * -1.0) * 10^6
  p7 <- ((x[11] > 2) * -1.0) * 10^6
  p16 <- ((x[16] > 2) * -1.0) * 10^6
  p17 <- ((x[17] > 2) * -1.0) * 10^6
  p18 <- ((x[18] > 2) * -1.0) * 10^6
  p1 <- ((x[1] > 10) * -1.0) * 10^6
  p17n <- ((x[17] < .0025) * -1.0) * 10^6
  prmin <- ((x[19] < 89014.12) * -1.0) * 10^6 #Del ultimo fit
  prmax <- ((x[19] > (263116 * 1.2)) * -1.0) * 10^6 #Del ultimo fit + 20%
  ptime <- (((1/x[5] + 1/x[6]) > 15) * -1.0) * 10^6 #Tiempo incubación/presintomático
  p5b <- ((x[9] < .1) * -1.0) * 10^6 #Tasa recuperación asintomáticos
  rm(z)
  rm(y)
  
  res + sum(penalty2 + penalty3 + penalty4 + penalty7 + penalty12 +
              penalty13 + penalty14 + penalty15 +
              penaltyb + p2 + p3 + p4 + p5 + p6 + p7 +
              p12 + p13 + p14 + p15 + p16 + p17 + p18 + p1 + p17n +
              prmin + prmax + ptime + p5b)
}


# Carga datos -------------------------------------------------------------

load("datos_CAN.RData")

f_inicio <- as.Date("2021-05-14")
f_fin <- as.Date("2021-06-12")

n_inicio <- which(data$FECHA == f_inicio)
n_fin <- which(data$FECHA == f_fin)

# Obtenemos datos diarios de fallecidos y número de tests
data$exitus <- c(0, diff(data$FALLECIDOS)) 
data$tests <- c(0, diff(data$TOTAL.TEST))

input2 <- data[n_inicio:n_fin, ]


# Estima desviación estandar de cada serie --------------------------------

# Esto nos permite asignar una likelihood al conjunto del modelo, considerando
# las series independientes.

devs <- c()
spl <- (input2$exitus)
fit <- smooth.spline(x = 1:nrow(input2), y = spl, df = 4)
devs[1] <- sd(spl - predict(fit)$y)
# plot(spl, type = "l")
# lines(predict(fit))

spl <- input2$TOTAL.HOSPITALIZADOS
fit <- smooth.spline(x = 1:nrow(input2), y = spl, df = 4)
devs[2] <- sd(spl - predict(fit)$y)

spl <- input2$Uci
fit <- smooth.spline(x = 1:nrow(input2), y = spl, df = 4)
devs[3] <- sd(spl - predict(fit)$y)

spl <- input2$HOSPITALIZADOS.UCI
fit <- smooth.spline(x = 1:nrow(input2), y = spl, df = 4)
devs[4] <- sd(spl - predict(fit)$y)

spl <- input2$Hosp
fit <- smooth.spline(x = 1:nrow(input2), y = spl, df = 4)
devs[5] <- sd(spl - predict(fit)$y)

spl <- (input2$CASOS.NUEVOS.PCR.)
fit <- smooth.spline(x = 1:nrow(input2), y = spl, df = 4)
devs[6] <- sd(spl - predict(fit)$y)


# Forzamientos ------------------------------------------------------------

forcings <- cbind(0:30, data$tests[n_inicio:(n_fin + 1)])

load("apple_CAN.RData")


aind1 <- cbind(0:30,
               apple.ind$unname.apple.[which(apple.ind$names.apple. == f_inicio):(which(apple.ind$names.apple. == f_fin) + 1)]/100)

forcings <- list(forcings, aind1)


# Hospitalizados iniciales ------------------------------------------------

hosp <- c()
hosp[1] <- data$TOTAL.HOSPITALIZADOS[n_inicio - 1]
hosp[2] <- data$HOSPITALIZADOS.UCI[n_inicio - 1]

# A punto de empezar el proceso -------------------------------------------

best <- -999999999 #LL inicial a mejorar

load("seeds_CAN.RData") #Cargamos seeds de los valores iniciales de los pars.

sols <- NA #Pre-aloco el número de combinaciones paramétricas en 2 unidades de LL de la mejor

set.seed(476468713)

condition <- T 
#Esta condición se hará falsa en el bucle de después, si en las rondas la LL no
#mejora más de 2 puntos y hay más de 1000 combinaciones paramétricas a menos de
#2 puntos de loglikelihood de la mejor solución.

round <- 1

# Ahora cambio los nombres de las columnas para que concuerden con los nombres
# en la función de likelihood
input2 <- input2 %>% rename(positius = CASOS.NUEVOS.PCR., 
                            ingressats = TOTAL.HOSPITALIZADOS,
                            novas_uci = Uci,
                            critics = HOSPITALIZADOS.UCI,
                            ingressos = Hosp)

# Paralelización

sims <- ncol(seeds) #Número de combinaciones paramétricas a explorar

Cores <- 19 #Numero de cores a utilizar.

while(condition){
  #Ahora viene la paralelización
  parall <- mclapply(1:sims, mc.cores = Cores, mc.preschedule = F,function(k){
    
    
    fit <- optim(par = seeds[, k], fn = ll_ode, forcings = forcings, y = input2, hosp = hosp, 
                 devs = devs, control = list(fnscale = -1, maxit = 500, parscale = seeds[, k]))
    
    if((k %% 1000) == 0) {
      cat("This was seed no. ", k, "\n")
      cat("This fit: ", fit$value, "\n")
    }
    
    fit
  })
  
  lhs <- parall
  
  rm(parall) #Para evitar fugas de memoria
  
  filename <- paste0("CAN_covid_20210612_", round, ".RData") #Salva cada ronda de optimizaciones, por si acaso
  save(lhs, file = filename)
  
  # Ahora, recuperamos la loglikelihood de cada combinación de parámetros
  logl <- rep(NA, 100000)
  for(i in 1:100000) logl[i] <- lhs[[i]]$value
  
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
  parmat <- matrix(NA, nrow = length(index), ncol = 19)
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
