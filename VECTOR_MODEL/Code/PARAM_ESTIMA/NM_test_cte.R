rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")

# Create Pseudo Data:
Path = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_vec_cte.c")
dyn.load("model_vec_cte.so")

f = 200
K = 250000
H = 160000
omega_t = 4
trueSD = 100
# We create a vector with the constant parameters.
parms = c(f,K,H,omega_t)
# We set the initial conditions to zero.
Y <- c(y1 = 100.0, y2 = 0.0, y3 = 0.0)
# List with the data frames of the forcings, sort as the c code.
min_t <- 1
max_t <- 365
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model_vec_cte",
           initfunc = "initmod", nout = 1,
           outnames = "Sum") 

ode <- data.frame(out) 
ode$Sum <- NULL

df_plot <- reshape2::melt(ode, id.vars = c("time"))
ggplot(df_plot,aes(time, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Counts") +
  ggtitle("Vector dynamics")+
  scale_color_manual(name = "",
                     labels = c("Larva", "Adult mosquito", "Adult handling mosquito"),
                     values=c('#FF00F6','#FF2C00','#FF2C23'))+
  theme_bw()

ggplot(ode) + 
  geom_line(aes(time,y3))

colnames(ode) <- c("time","L","A","Ah")
head(ode)

omega_t = 1
trueSD = 100
# We create a vector with the constant parameters.
parms = c(f,K,H,omega_t)
# We set the initial conditions to zero.
Y <- c(y1 = 100.0, y2 = 0.0, y3 = 0.0)
times <- seq(min_t,max_t, 1)
out2 <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model_vec_cte",
           initfunc = "initmod", nout = 1,
           outnames = "Sum") 

ode2 <- data.frame(out2) 
ode2$Sum <- NULL

df_plot2 <- reshape2::melt(ode2, id.vars = c("time"))
ggplot(df_plot2,aes(time, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Counts") +
  ggtitle("Vector dynamics")+
  scale_color_manual(name = "",
                     labels = c("Larva", "Adult mosquito", "Adult handling mosquito"),
                     values=c('#FF00F6','#FF2C00','#FF2C23'))+
  theme_bw()

ggplot(ode2) + 
  geom_line(aes(time,y3))

colnames(ode2) <- c("time","L","A","Ah")
sum1 <- ode$A + ode$Ah
sum2 <- ode2$A + ode$Ah
head(ode2)
diff_df <- abs(sum1 - sum2)
saveRDS(ode, file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo.rds")
# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE
Path = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model_vec_cte.c")
dyn.load("model_vec_cte.so")


# Función que calcula la loglikelihood del modelo ----------------------------------

# Esta será la función a optimizar

ll_ode <- function(x, # Params
                   y, # datos
                   devs){ #desviaciones estándar para calcular la loglikelihood
  
  # if(x[1] <=  0){
  if(1 ==  0){
    res = -86829146000
  }else{
    pars <- c(f = f,K = K,H = H,omega = x[1])
    cat("Pars:",pars, "\n")
    population <- c(y1 = 100.0, y2 = 0.0, y3 = 0.0)
    
    z <- ode(y=population,
             times = 0:nrow(y), func = "derivs", 
             dllname = "model_vec_cte" , parms = pars,
             initfunc = "initmod") 
    
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
  cat("res:",res, "\n")
  return(res)
}

# Carga datos -------------------------------------------------------------
# Pseudo Data to check the oprimization method.
ob_data <- readRDS(file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo.rds")
colnames(ob_data) <- c("time", "L", "A", "Ah")
l <- nrow(ob_data)
ob_data$L <- ob_data$L + rnorm(l,0,trueSD)
ob_data$A <- ob_data$A + rnorm(l,0,trueSD)
ob_data$Ah <- ob_data$Ah + rnorm(l,0,trueSD)
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

y <- ode
# (x, # vector con los parámetros
#   forcs_mat, # forzamientos para el solver de la ode
#   y, # datos
#   devs)
likelyhood <- function(x){
  return(ll_ode(x,y,devs))
}
x <- seq(1, 500, by=.05)
slopelikelihoods <- lapply(x,likelyhood)
slopelikelihoods_num <- as.numeric(unlist(slopelikelihoods))
slopelikelihoods_num <- slopelikelihoods_num - trunc(slopelikelihoods_num)

plot (x, slopelikelihoods , type="l", xlab = "values of omega", ylab = "Log likelihood")# Forzamientos ------------------------------------------------------------
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
