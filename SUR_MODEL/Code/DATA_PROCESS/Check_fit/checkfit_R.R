rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")

# Read output from param_estimation.R
output <- load("/home/marta/Documentos/PHD/2021/SUR_Model/OUTPUT/param_MAD_MODEL_1core_900it1.RData")
# Number of seeds:
N = 300
out_mat <- matrix(0, nrow = N, ncol = 4)
for(i in c(1:N)){
  out_mat[i,1] = lhs[[i]]$value
  out_mat[i,2:4] = lhs[[i]]$par
}

mat_sort <- out_mat[order(out_mat[,1], decreasing = T),]

# Load C code:
Path = "~/MAD_MODEL/SUR_MODEL/Code/PARAM_ESTIMATION/MH/"
setwd(Path)
system("R CMD SHLIB model_2000eq.c")
dyn.load("model_2000eq.so")

# Real value parameters:
it <- 1
gam1 = mat_sort[it,2]
gam2 = mat_sort[it,3]
gam3 = mat_sort[it,4]

pars <- c(gam1,gam2,gam3)
# Read Observed data:
ob_data <- read.table("~/MAD_MODEL/SUR_MODEL/data/Observed_2021-11-23_2691.dat", header=FALSE, sep= "\t")
ob_data_t <- t(ob_data)
downloads <- ob_data_t[,1:2]
end_time <- max(ob_data_t[,1])
age_max <- max(ob_data_t[1,]) -1
population <- matrix(0,nrow = 1,ncol=2000)
z <- ode(y = population,
         times = 0:end_time, func = "derivs", method = "ode45",
         dllname = "model_2000eq", initfunc = "initmod", nout = 0, 
         parms = pars, initforc = "forcc", forcings = downloads, 
         fcontrol = list(method = "constant"))

z <- as.data.frame(z)
