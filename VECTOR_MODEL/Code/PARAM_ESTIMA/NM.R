rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")

Path_temp = "~/MAD_MODEL/VECTOR_MODEL/data/bcn_weather_daily.Rds"
# Path_temp = paste(PC,Path_temp, sep="")

temp <-read_rds(Path_temp)
temp$date = as.Date(temp$date , "%Y-%m-%d")
temp <- temp %>%  group_by(date) %>% summarise(mean_temp = mean(valor))
# 
# ggplot(temp) +
#   geom_line(aes(date, mean_temp))+
#   ggtitle("Mean temperature Barcelona") + 
#   xlab("Mean temperature")+
#   theme_bw()

# Gonotrophic cycle:
gonot <- function(T){
  val = (0.045*T^2 - 2.717*T + 44.405)^(-1)
  return(val)
}

# Development rate
d_L <- function(T){
  sigma <- (0.14457*T^2 - 8.24857*T + 124.80857)^(-1)
  return(sigma)
}

# Larva mortality rate
delta_L <- function(T){
  rate <- 0.021643*T^2 - 0.959568*T + 10.440131 
  if (rate < 0) {
    rate <- 0 }
  return(rate)
}

# Adult mosquito mortality rate
delta_A <- function(T){
  if (T <=20) {
    rate <- - 0.01019105*T^3 + 0.5920223*T^2 - 11.38171*T + 72.60858
    #rate <- 0
  } else if (T>20 & T<=25) {
    rate <- 0.001216109*T^3 - 0.09240726*T^2 + 2.306881*T - 18.6487 #rate <- 0
  } else if (T>25 & T<=30) {
    rate <- 0.002551257*T^3 - 0.1925433*T^2 + 4.810282*T - 39.51037
  } else
    rate <- 0.01310967*T^3 - 1.142801*T^2 + 33.31801*T - 324.5876
  return(rate)
}



vec = seq(0,40,1)
df_gonot_vec <- data.frame(temp = vec, gonot =unlist(lapply(vec,gonot)))
df_dL_vec <- data.frame(temp = vec, dL = unlist(lapply(vec,d_L)))
df_deltaL_vec <- data.frame(temp = vec, deltaL = unlist(lapply(vec,delta_L)))
df_deltaA_vec <- data.frame(temp = vec, deltaA = unlist(lapply(vec,delta_A)))
# df_dl_opt_vec <- data.frame(temp = vec, deltaA = unlist(lapply(vec,dL_opt)))
# 
# ggplot(df_gonot_vec) + geom_line(aes(temp,gonot)) + 
#   ggtitle("Inverse of the Gonotrophic cycle")+
#   theme_bw()
# 
# ggplot(df_dL_vec) + geom_line(aes(temp,dL)) +
#   ggtitle("Larva development rate")+
#   theme_bw()
# 
# ggplot(df_deltaL_vec) + geom_line(aes(temp,deltaL)) +
#   ggtitle("Larva mortality rate")+
#   theme_bw()
# 
# ggplot(df_deltaA_vec) + geom_line(aes(temp,deltaA)) +
#   ggtitle("Mosquito adult mortality rate")+
#   theme_bw()

# ggplot(df_dl_opt_vec) + geom_line(aes(temp,deltaA)) +
#   ggtitle("Mosquito adult mortality rate")+
#   theme_bw()


# Compute the values of the functions/forcings with temp.
# Compute the minimum date of the rho:
min_date <- min(temp$date)
max_date <- max(temp$date)
# DFs with the date and value of the parameter at that time.
temp <- temp  %>% filter( temp$date >= min_date & temp$date <= max_date)
df_date <- data.frame(date = temp$date)
df_date$time = as.numeric(df_date$date - as.Date(min_date,"%Y-%m-%d") , units="days") 

gono = unlist(lapply(temp$mean_temp,gonot))
df_gonot_out <- data.frame(date = temp$date, gono)

df_gonot_out$time = as.numeric(df_gonot_out$date - as.Date(min_date,"%Y-%m-%d") , units="days") 
df_gonot_out <- df_gonot_out %>% filter( df_gonot_out$time >= 0)

df_dL_out <- data.frame(date = temp$date, dL = unlist(lapply(temp$mean_temp,d_L)))
df_dL_out$time = as.numeric(df_dL_out$date - as.Date(min_date,"%Y-%m-%d") , units="days") 
df_dL_out <- df_dL_out %>% filter( df_dL_out$time >= 0)

df_deltaL_out <- data.frame(date = temp$date, deltaL = unlist(lapply(temp$mean_temp,delta_L)))
df_deltaL_out$time = as.numeric(df_deltaL_out$date - as.Date(min_date,"%Y-%m-%d") , units="days") 
df_deltaL_out <- df_deltaL_out %>% filter( df_deltaL_out$time >= 0)

df_deltaA_out <- data.frame(date = temp$date, deltaA = unlist(lapply(temp$mean_temp,delta_A)))
df_deltaA_out$time = as.numeric(df_deltaA_out$date - as.Date(min_date,"%Y-%m-%d") , units="days") 
df_deltaA_out <- df_deltaA_out %>% filter( df_deltaA_out$time >= 0)
# 
# ggplot(df_gonot_out) +
#   geom_line(aes(date,gono)) +
#   ggtitle("Gonotrophic cycle")+
#   theme_bw()
# 
# ggplot(df_dL_out) +
#   geom_line(aes(date,dL)) +
#   ggtitle("Larva development rate")+
#   theme_bw()
# 
# ggplot(df_deltaL_out) +
#   geom_line(aes(date,deltaL)) +
#   ggtitle("Larva mortality rate")+
#   theme_bw()
# # 
# ggplot(df_deltaA_out) +
#   geom_line(aes(date,deltaA)) +
#   ggtitle("Adult mosquito mortality rate")+
#   theme_bw()

df_gonot_out$date <- NULL
df_gonot_out <- df_gonot_out[,c(2,1)]
head(df_gonot_out)
df_dL_out$date <- NULL
df_dL_out <- df_dL_out[,c(2,1)]
head(df_dL_out)
df_deltaL_out$date <- NULL
df_deltaL_out <- df_deltaL_out[,c(2,1)]
head(df_deltaL_out)
df_deltaA_out$date <- NULL
df_deltaA_out <- df_deltaA_out[,c(2,1)]
head(df_deltaA_out)


# Create pseudo data:
# Create Pseudo Data:
Path = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model.c")
dyn.load("model.so")

f = 200
K = 250000
H = 160000
omega_t = 4
trueSD = 1
# We create a vector with the constant parameters.
parms = c(f,K,H, omega_t)
# We set the initial conditions to zero.
Y <- c(y1 = 100.0, y2 = 0.0, y3 = 0.0)
# List with the data frames of the forcings, sort as the c code.
forcs_mat <- list(data.matrix(df_gonot_out),
                  data.matrix(df_dL_out),
                  data.matrix(df_deltaL_out),
                  data.matrix(df_deltaA_out))
min_t <- min(df_dL_out$time)
max_t <- max(df_dL_out$time)
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model",
           initfunc = "initmod", nout = 1,
           outnames = "Sum", initforc = "forcc",
           forcings = forcs_mat, fcontrol = list(method = "constant")) 

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

saveRDS(ode, file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo.rds")
# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE
Path = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model.c")
dyn.load("model.so")


# Función que calcula la loglikelihood del modelo ----------------------------------

# Esta será la función a optimizar

ll_ode <- function(x, # vector con los parámetros
                   forcs_mat, # forzamientos para el solver de la ode
                   y, # datos
                   devs){ #desviaciones estándar para calcular la loglikelihood
  
  if(x[1] <=  0){
    res = -86829146000
  }else{
    pars <- c(f = f,K = K,H = H,omega = x[1])
    
    population <- c(y1 = 100.0, y2 = 0.0, y3 = 0.0)
    
    z <- ode(y=population,
             times = 0:nrow(y), func = "derivs", 
             dllname = "model" , parms = pars,
             initfunc = "initmod", initforc = "forcc",
             forcings = forcs_mat, fcontrol = list(method = "constant")) 
    #Aquí corre el ODE
    
    colnames(z)[2:4] <- c("L", "A", "Ah")
    
    z <- as.data.frame(z)
    z <- z[-1, ]
    
    L <- y$L
    A <- y$A
    Ah <- y$Ah
    
    res <- #cálculo de la loglikelihood en función de las desviaciones estándar
      sum(dnorm((A+Ah), mean = (z$A + z$Ah), sd = devs[1], log = T))  
  }
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
# head(ob_data)

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


# Forzamientos ------------------------------------------------------------

# Registrations:
# head(down)

# Hospitalizados iniciales ------------------------------------------------

# A punto de empezar el proceso -------------------------------------------

best <- -999999999 #LL inicial a mejorar

# load("seeds_CAN.RData") #Cargamos seeds de los valores iniciales de los pars.
seeds <- matrix(runif(10000,0,1), ncol = 10000, nrow = 1)
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

Cores <- parallel::detectCores() #Numero de cores a utilizar.
it <- 0
while(condition){
  #Ahora viene la paralelización
  parall <- mclapply(1:sims, mc.cores = Cores, mc.preschedule = F,function(k){
    it <- it + 1
    
    fit <- optim(par = seeds[, k], fn = ll_ode, forcs_mat = forcs_mat, y = input2, 
                 devs = devs, control = list(fnscale = -1, maxit = 500, parscale = seeds[, k]))
    
    if((k %% 1000) == 0) {
      cat("This was seed no. ", k, "\n")
      cat("This fit: ", fit$value, "\n")
    }
    
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
