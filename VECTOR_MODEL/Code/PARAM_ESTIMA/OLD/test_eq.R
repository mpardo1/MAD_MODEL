rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("coda")

# Params values:
a = 2 
b = 3
c = 2.8
d = 5
# Equilibrium points:
x_eq <- -c/a
y_eq <- 0

###############   ODE INTEGRATION   ##################
Path = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB test_mod.c")
dyn.load("test_mod.so")


trueSD = 1
# We create a vector with the constant parameters.
parms = c(a,b,c,d)
# We set the initial conditions to zero.
Y <- c(y1 = x_eq, y2 = y_eq)
min_t <- 1
max_t <- 4200
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs", method = "ode45",
           parms = parms, dllname = "test_mod",
           initfunc = "initmod", nout = 1,
           outnames = "Sum") 

ode <- data.frame(out) 
ode$Sum <- NULL


# saveRDS(ode, file = "~/MAD_MODEL/VECTOR_MODEL/Code/PARAM_ESTIMA/ode_pseudo_cte.rds")
df_plot <- reshape2::melt(ode, id.vars = c("time"))
ggplot(df_plot,aes(time, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Counts") +
  ggtitle("Vector dynamics")+
  scale_color_manual(name = "",
                     labels = c("Larva", "Adult mosquito", "Adult handling mosquito"),
                     values=c('#FF00F6','#FF2C00','#2F822B'))+
  theme_bw()

