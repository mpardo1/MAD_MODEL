rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("coda")
#-----------------------------------------------------------------------------------#
# Equilibrium points:
eq_point <- function(delta_L, delta_A, d_L, a, fec, K, Hum, omega_t){
  Ah_eq <-  K*(((omega_t*Hum*d_L)/((a+delta_A)*(omega_t*Hum +delta_A)))-((1/(a*fec))*(d_L+delta_L)))
  L_eq <- ((a*fec*Ah_eq)/((a*fec/K)*Ah_eq+(d_L+delta_L)))
  A_eq <- ((d_L/(omega_t*Hum+delta_A))*((a*fec*Ah_eq)/(((a*fec/K)*Ah_eq)+(d_L+delta_L))))
  eq <- c(L_eq, A_eq, Ah_eq)
}

# Stability condition:
feasability_cond <- function(del_L, del_A, dev_L, gon, fec, Kar, Hum, omega){
  fes <- FALSE
  val <- (gon*fec*omega*Hum*dev_L)/((gon + del_A)*(omega*Hum + del_A)*(dev_L + del_L))
  if(val > 1){
    fes <- TRUE
  }
  return(fes)
}

# ODE system in R:
vect <- function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    
    dL <-  gon*fecun*H*(1-(L/Ka))-(dev_L+del_L)*L
    dA <-  dev_L*L - (omeg + del_A)*A
    dH <-  omeg*A - (gon + del_A)*H
    # return the rate of change
    list(c(dL, dA, dH))
  }) # end with(as.list ...
}


# Params values:
fec = 2000
K = 4300
Hum = 13554
omega_t = 0.2
delta_L = 0.2
delta_A = 0.3
d_L = 14
a = 0.01

vec_eq <- eq_point(delta_L, delta_A, d_L, a, fec, K, Hum, omega_t)

times <- seq(0, 15, by = 0.01)
parameters <- c(fecun = fec,Ka = K,Hu = Hum,omeg = omega_t,del_L = delta_L,del_A = delta_A,dev_L = d_L,gon = a)
# state <- c(L = vec_eq[1], A = vec_eq[2], H = vec_eq[3])
state <- c(L = 10, A = 0, H = 0)
out2 <- as.data.frame(ode(y = state, times = times, func = vect, parms = parameters))
# df_plot2 <- reshape2::melt(out2, id.vars = c("time"))
# ggplot(df_plot2,aes(time, value))  +
#   geom_line(aes( colour = variable)) +
#   ylab("Counts") +
#   ggtitle("Vector dynamics")+
#   scale_color_manual(name = "",
#                      labels = c("Larva", "Adult mosquito", "Adult handling mosquito"),
#                      values=c('#FF00F6','#FF2C00','#2F822B'))+
#   theme_bw()
saveRDS(out2, file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo_cte.rds")
# Función que calcula la loglikelihood del modelo ----------------------------------

# Esta será la función a optimizar

ll_ode <- function(x, # Params
                   y, # datos
                   devs){ #desviaciones estándar para calcular la loglikelihood
  
 if(x <=  0){
    res = -86829146000
  }else{
    pars <- c(fecun = fec,Ka = K,Hu = Hum,omeg = x,
              del_L = delta_L,del_A = delta_A,dev_L = d_L,
              gon = a) # death rate group 2
    
    cat("Pars:",pars, "\n")
    population <- c( 10,0,0)
    
    z <-  as.data.frame(ode(y = state, times = times, func = vect, parms = parameters))
    
    #Aquí corre el ODE
    
    colnames(z)[2:4] <- c("L", "A", "Ah")
    
    z <- as.data.frame(z)
    z <- z[-1, ]
    
    L <- y$L
    A <- y$A
    Ah <- y$Ah
    
    res <- sum(dnorm(L, mean = z$L, sd = 1, log = T))  +
           sum(dnorm(A, mean = z$A , sd = 1, log = T))  +
           sum(dnorm(Ah, mean = z$Ah, sd = 1, log = T))  
  }
  # cat("res:",res, "\n")
  return(res)
}

# Carga datos -------------------------------------------------------------
# Pseudo Data to check the oprimization method.
ob_data <- readRDS(file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo_cte.rds")
colnames(ob_data) <- c("time", "L", "A", "Ah")
l <- nrow(ob_data)
# ob_data$L <- ob_data$L + rnorm(l,0,trueSD)
# ob_data$A <- ob_data$A + rnorm(l,0,trueSD)
# ob_data$Ah <- ob_data$Ah + rnorm(l,0,trueSD)
head(ob_data)
df_plot <- reshape2::melt(ob_data, id.vars = c("time"))
ggplot(df_plot,aes(time, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Counts") +
  ggtitle("Vector dynamics")+
  scale_color_manual(name = "",
                     labels = c("Larva", "Adult mosquito", "Adult handling mosquito"),
                     values=c('#FF00F6','#FF2C00','2F822B'))+
  theme_bw()



# f_inicio <- as.Date("2021-05-14")
# f_fin <- as.Date("2021-06-12")
# 
# n_inicio <- which(data$FECHA == f_inicio)
# n_fin <- which(data$FECHA == f_fin)

n_inicio <- 1
n_fin <- max(ob_data$time)

input2 <- ob_data[n_inicio:n_fin, ]


# Estima desviación estandar de cada serie --------------------------------

# Esto nos permite asignar una likelihood al conjunto del modelo, considerando
# las series independientes.

devs <- c()
spl <- input2$A + input2$Ah
fit <- smooth.spline(x = 1:nrow(input2), y = spl, df = 4)
devs[1] <- sd(spl - predict(fit)$y)

y <- as.data.frame(ob_data)
# (x, # vector con los parámetros
#   forcs_mat, # forzamientos para el solver de la ode
#   y, # datos
#   devs)
likelyhood <- function(x){
    return(ll_ode(x,y,devs))
  }

x <- seq(1.e-5, 1, by=.0001)
slopelikelihoods <- lapply(x,likelyhood)
slopelikelihoods_num <- as.numeric(unlist(slopelikelihoods))
slopelikelihoods_num <- slopelikelihoods_num - trunc(slopelikelihoods_num)

plot (x, slopelikelihoods_num , type="l", xlab = "values of omega", ylab = "Log likelihood")# Forzamientos ------------------------------------------------------------
plot (x, slopelikelihoods , type="l", xlab = "values of omega", ylab = "Log likelihood")# Forzamientos ------------------------------------------------------------

# Registrations:
# head(down)

# Hospitalizados iniciales ------------------------------------------------

# A punto de empezar el proceso -------------------------------------------

best <- -999999999 #LL inicial a mejorar
# We set the seed if we want to replicate the results again.
#set.seed(4468713)
# load("seeds_CAN.RData") #Cargamos seeds de los valores iniciales de los pars.
seeds <- matrix(runif(10000,0,4), ncol = 10000, nrow = 1)
sols <- NA #Pre-aloco el número de combinaciones paramétricas en 2 unidades de LL de la mejor



condition <- T 
#Esta condición se hará falsa en el bucle de después, si en las rondas la LL no
#mejora más de 2 puntos y hay más de 1000 combinaciones paramétricas a menos de
#2 puntos de loglikelihood de la mejor solución.

round <- 1
# Paralelización

# sims <- 2 #Número de combinaciones paramétricas a explorar
sims <- ncol(seeds) #Número de combinaciones paramétricas a explorar

Cores <- 16 #Numero de cores a utilizar.
it <- 0
while(condition){
  #Ahora viene la paralelización
  parall <- lapply(1:sims,function(k){
    it <- it + 1
    
    fit <- optim(par = seeds[, k], fn = ll_ode, y = input2, 
                 devs = devs, control = list(fnscale = -1, maxit = 500, parscale = seeds[, k]))
    
    if((k %% 1000) == 0) {
      cat("This was seed no. ", k, "\n")
      cat("This fit: ", fit$value, "\n")
    }
    cat("This was seed no. ", k, "\n")
    cat("This fit: ", fit$value, "\n")
    cat("The proposal is:",fit$par ,"\n")
    fit
  })
  
  lhs <- parall
  
  rm(parall) #Para evitar fugas de memoria
  
  filename <- paste0("~/MAD_MODEL/VECTOR_MODEL/OUTPUT/param_VECTOR_MODEL", round, ".RData") #Salva cada ronda de optimizaciones, por si acaso
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
