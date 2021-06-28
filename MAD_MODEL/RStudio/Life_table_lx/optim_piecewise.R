rm(list = ls())
library(gdata) 
library(segmented)
library(e1071)
library(ggplot2)
library(optimx)
library(zoom)

# Read data from de csv and save it on life_table:
life_table = read.csv("Documentos/R/data/participation_life_table.csv") 

# Segmented does not allow referenced to the data frame class
y <- log10(life_table$lx/life_table$lx[1])
x <- life_table$x 

# With Unkown thresholds. 
gamma_func <- function(data, par){
  with(data, sum((  (x<=par[1])*(log10( par[2]/(par[3] + par[2])) *(x - 1) )
                  + (x>par[1] & x<=par[4])*(log10( par[5]/(par[6] + par[5])) *(x - 1) )
                  #+ (x>par[4] & x<=par[7])*(log10( par[8]/(par[8] + par[9])) *(x - 1) ) 
                  - y)^2))
}

dat = data.frame(x=life_table$x,y=log10(life_table$lx/life_table$lx[1]))
result <- optim(par = c(1,2,3,4,5,6,7,8,9), fn = gamma_func, lower= c(0,1.e-6,1.e-6,0,1.e-6,1.e-6,0,1.e-6,1.e-6), method="L-BFGS-B", data = dat)
result_1 <- optim(par = c(1,2,3,4,5,6), fn = gamma_func, lower= c(0,1.e-6,1.e-6,0,1.e-6,1.e-6), method="L-BFGS-B", data = dat)


param = result$par
x <- dat$x
c1 = param[2]/(param[3] + param[2])
c2 = param[5]/(param[5] + param[6])
c3 = param[8]/(param[8] + param[9])
func = (x<param[1])*(log10( c1) *(x - 1) ) + (x>param[1] & x<=param[4])*(log10( c2) *(x - 1) ) + (x>param[4] & x<=param[7])*(log10( c3) *(x - 1) ) 
func = (x<param[1])*(log10( c1) *(x - 1) ) + (x>param[1] & x<=param[4])*(log10( c2) *(x - 1) ) 

plot(dat$x,dat$y, main="Least square regression")
lines(x,func)


#With known thresholds
gamma_func1 <- function(data, par){
  with(data, sum((  (x<=105)*(log10( par[1]/(par[1] + par[2])) *(x - 1) )
                    + (x>'105')*(log10( par[3]/(par[3] + par[4])) *(x - 1) )
                    #+ (x>par[4] & x<=par[7])*(log10( par[8]/(par[8] + par[9])) *(x - 1) ) 
                    - y)^2))
}

dat1 = data.frame(x=life_table$x,y=log10(life_table$lx/life_table$lx[1]))
result1 <- optim(par = c(1,2,3,4), fn = gamma_func1, lower= c(0,0,0,0), method="L-BFGS-B", data = dat1)


param1 = result1$par
x1 <- dat$x
c11 = param1[1]/(param1[1] + param1[2])
c21 = param1[3]/(param1[3] + param1[4])
func1 = (x<=105)*(log10( c11) *(x - 1) ) + (x>105)*(log10( c21) *(x - 1) ) 

plot(dat$x,dat$y, main="Least square regression")
lines(x1,func1)
