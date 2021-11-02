rm(list = ls())
library(gdata) 
library(segmented)
library(e1071)
library(ggplot2)
library(optimx)
library(zoom)
library(deSolve)
library(mlr3misc)
library(numbers)

# Read data from de csv and save it on life_table:
download = read.csv("users.csv") 
download$date_trunc = strtrim(download$date, 10)
download$date_1 = as.Date(download$date_trunc, "%Y-%m-%d")
download$time_down = seq(from = 1, to = 7*length(download$date), by = 7)
download$mosquito_count = 1

for( i in c(2:length(download$date))){
  download$mosquito_count[i] = download$n[i] - download$n[i-1]
}
lenght_days = length(download$time_down)
downloads_uni = data.frame(days = c(1:download$time_down[lenght_days]), number_down = c(1:download$time_down[lenght_days]))
downloads_uni$number_down[1:7] = download$mosquito_count[1]/7
length_data = length(download$n)-2
for(i in c(1:length_data)){
  init = 7*i
  out = init + 7
  # browser() Para debuggear.
  downloads_uni$number_down[init:out] = download$mosquito_count[i+1]/7
 # downloads_uni$number_down[7*i:7*i+7] = download$mosquito_count[i+1]/7
}
downloads_uni$number_down[length(downloads_uni$days)] = download$mosquito_count[length_data+2]/7
# Start with a toy model coded in R
# pars_simp <- list(A = 200, mu = 0.5, gamma = 0.3)
# 
# solveP_simp <- function(pars_simp, times=seq(0,50,by=0.5)) {
#   derivs <- function(t, state, pars_simp) { # returns rate of change
#         with(as.list(c(state, pars_simp)), {
#             dP1 <- A - mu*P1 - gamma*P1
#             dP2 <- mu*P1 - mu*P2 - gamma*P2
#             dP3 <- mu*P2 - mu*P3
#                 return(list(c(dP1,dP2,dP3), TOC = P1 + P2 + P3))
#         })
#     }
#       state<- c(P1 = 100, P2 = 100, P3 = 100) ## ode solves the model by integration...
#         return(ode(y = state, times = times, func = derivs, parms = pars_simp))
# }
# 
# out_simp <- solveP_simp(pars_simp)
# 
# matplot(out_simp[,1], out_simp[,-1], type = "l", lty = 1:3, lwd = c(2, 2, 1), col = "black", xlab = "time, hour", ylab = "Total number")
# legend("topright", c("P1", "P2", "P3", "TOC"),lty = 1:3, lwd = c(2, 2, 2, 1))


# Number of ages + 1 (Number of ages should match with C variable = NUM_AGE)
num_age = 2000+1
# Make sure you are in the working directory of the file toy_model.c
# If not uncomment the following two lines setting Path = the path where the toy_model.c file is saved.
Path = "~/Documentos/SUR Model/Code/RStudio/Fitting"
setwd(Path)
system("R CMD SHLIB toy_model.c")
dyn.load("toy_model.so")
# Toy model code in C.
forcs_mat = data.matrix(downloads_uni)
# This should coincide with the parameters specified in the C code.
parms <- c(mu1  = 0.04, gamma1 = 0.4, mu2= 0 ,gamma2=0, mu3 = 0, gamma3 = 0)
Y <- numeric(num_age)
#Y  <- c(y1 = 0, y2 = 0.0, y3 = 0.0, y4 = 0, y5 = 0)
times <- 0:100
out <- ode(Y, times, func = "derivs", parms = parms,
           dllname = "toy_model", initfunc = "initmod", nout = 1, outnames = "Sum",
           initforc = "forcc", forcings = forcs_mat, fcontrol = list(method = "constant")
)
end_vec = length(out[1,])-1
start_vec = length(out[1,])-3
out[,start_vec:end_vec]
# Prove visually that the forcing is working, i.e., the peaks on the first group ages are at the peaks of downloads
out_df <- data.frame( Time = out[,1], P1 = out[,2], P2 = out[,3])
ggplot( out_df, aes(x = Time, y = P1)) + geom_line() + 
  geom_point(data = downloads_uni, aes(x = days, y = number_down), color = "red")

# The download function is piecewise since it is the uniform distribution at each week, it can bee seen at shorter scale
plot(downloads_uni$days[1:10],downloads_uni$number_down[1:10])

