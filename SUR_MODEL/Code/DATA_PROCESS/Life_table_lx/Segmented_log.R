rm(list = ls())
library(gdata) 
library(segmented)
library(e1071)
library(ggplot2)


# Read data from de csv and save it on life_table:
life_table = read.csv("Documentos/R/data/participation_life_table.csv") 
lx <- log10(life_table$lx/(life_table$lx[1]))
x <- life_table$x
#First we create a linear model
LinearMod <- lm(lx~x)

dat_nls = data.frame(x=life_table$x,y=log10(life_table$lx/life_table$lx[1]))
LinearMod1 <- nls(lx~(log10(mu/(mu+gamma))*x-log10(mu/(mu+gamma))), data = dat_nls, start = list(mu = 0.1, gamma = 0.01))
#Create vectors to save the squared errors and the number of breakpoints.
rSq_ve = c(1:10)
nod_vec = c(1:10)

# Create a loop over the fixed number of breakpoints for the segmented algorithm.
for(i in 1:10){
  seg <- segmented(LinearMod, seg.Z = ~x, npsi = i)
  nod_vec[i] = i
  rSq_ve[i] = summary(seg)$adj.r.squared
}

#Plot the squared error against the number of breakpoints.
plot(nod_vec, rSq_ve, xlab = "Number of breakpoints", ylab = "R squared", xlim = c(0,10))

# Do segmented for 5 breakpoints since it appears to be the most efficient number of breakpoints.
num_breakp = 3
#seg_5 <- segmented(LinearMod, seg.Z = ~x, control=seg.control( fix.npsi=FALSE))
seg_5 <- segmented(LinearMod, seg.Z = ~x, npsi = num_breakp)
fit_data2 <- fitted(seg_5)
model <- data.frame(Age = x, Participants = fit_data2)
ggplot(model, aes(x = Age, y = Participants)) + geom_line() 
ggplot(model, aes(x = Age, y = Participants)) + geom_line() + geom_point(data = life_table, aes(x = x, y = log10(lx/lx[1])), color = "red")

gamma_func <- function(data, par){
    with(data, sum((log10( par[1]/(par[1] + par[2])) *x - log10( par[1]/(par[1] + par[2]))  - y)^2))
}
dat = data.frame(x=life_table$x[1:16],y=log10(life_table$lx[1:16]/life_table$lx[1]))

result <- optim(par = c(1,2), fn = gamma_func, lower= c(0,0), method="L-BFGS-B", data = dat)
result

plot(y ~ x, data = dat, main="Least square regression")
abline(a = result$par[1], b = result$par[2], col = "red")

dim_vec = 16
age_max = 107
x_vec = c(1:age_max)
p_vec = c(1:age_max)
mu = result$par[1]
gamma = result$par[2]
g = 0
for(i in 1:age_max){
  g = g + log10(mu/(mu + gamma))
  p_vec[i] = g
}

plot(x_vec,p_vec, pch = 18)
points(life_table$x[1:dim_vec], log10(life_table$lx[1:dim_vec]/life_table$lx[1]), col = 'green')



dat1 = data.frame(x=life_table$x[16:155],y=log10(life_table$lx[16:155]/life_table$lx[1]))

result1 <- optim(par = c(1,2), fn = gamma_func, lower= c(0,0), method="L-BFGS-B", data = dat1)
result1


age_min = 107
age_max = 1079
dim_vec = length(life_table$x[age_min:age_max])
x1 = c(age_min:age_max)
p_vec1 = c(age_min:age_max)
mu = result1$par[1]
gamma = result1$par[2]
g = 0
for(i in 1:age_max-age_min){
  g = g + log10(mu/(mu + gamma))
  p_vec1[i] = g
}

plot(x1,p_vec1, pch = 18)
points(life_table$x, log10(life_table$lx/life_table$lx[1]), col = 'green')

dat2 = data.frame(x=life_table$x[155:174],y=log10(life_table$lx[155:174]/life_table$lx[1]))

result2 <- optim(par = c(1,2), fn = gamma_func, lower= c(0,0), method="L-BFGS-B", data = dat2)
result2

age_min = 1079
age_max = 1216
dim_vec = age_max-age_min
x2 = c(age_min:age_max)
p_vec2 = c(age_min:age_max)
mu = result2$par[1]
gamma = result2$par[2]
g = 0
for(i in 1:age_max-age_min){
  g = g + log10(mu/(mu + gamma))
  p_vec2[i] = g
}

plot(x2,p_vec2, pch = 18)
points(life_table$x[155:174], log10(life_table$lx[155:174]/life_table$lx[1]), col = 'green')

dat3 = data.frame(x=life_table$x[174:254],y=log10(life_table$lx[174:254]/life_table$lx[1]))

result3 <- optim(par = c(1,2), fn = gamma_func, lower= c(0,0), method="L-BFGS-B", data = dat3)
result3

age_min = 1216 
age_max = 1771
dim_vec = age_max-age_min
x3 = c(age_min:age_max)
p_vec3 = c(age_min:age_max)
mu = result3$par[1]
gamma = result3$par[2]
g = 0
for(i in 1:age_max-age_min){
  g = g + log10(mu/(mu + gamma))
  p_vec3[i] = g
}

plot(x3,p_vec3, pch = 18)
points(life_table$x[174:254], log10(life_table$lx[174:254]/life_table$lx[1]), col = 'green')

dat4 = data.frame(x=life_table$x[254:283],y=log10(life_table$lx[254:283]/life_table$lx[1]))

result4 <-  optim(par = c(1,2), fn = gamma_func, lower= c(0.1,0), method="L-BFGS-B", data = dat4)
result4

age_min = 1771 
age_max = 1974
dim_vec = age_max-age_min
x4 = c(age_min:age_max)
p_vec4 = c(age_min:age_max)
mu = result4$par[1]
gamma = result4$par[2]
g = 0
for(i in 1:age_max-age_min){
  g = g + log10(mu/(mu + gamma))
  p_vec4[i] = g
}

plot(x4,p_vec4, pch = 18)
points(life_table$x[254:1974], log10(life_table$lx[254:1974]/life_table$lx[1]), col = 'green')

age_min = 1
age_max = 1216
x_tot = c(age_min:age_max)
p_vec_tot = c(age_min:age_max)
g = 0
for(i in 1:age_max){
  if( i < 107){
    mu = result$par[1]
    gamma = result$par[2]
    print("i",1)
  }else if( i >= 107 & i < 1079){
    mu = result1$par[1]
    gamma = result1$par[2]
    print("2")
  }else if(i >= 1079 & i < 1216 ){
    mu = result2$par[1]
    gamma = result2$par[2]
    print("3")
  }else  if(i >= 1216 & i < 1771){
    mu = result3$par[1]
    gamma = result3$par[2]
    print("4")
  }else{
    mu = result4$par[1]
    gamma = result4$par[2]
  }
  g = g + log10(mu/(mu + gamma))
  p_vec_tot[i] = g
}

plot(x_tot,p_vec_tot, pch = 18)
points(life_table$x, log10(life_table$lx/life_table$lx[1]), col = 'green')

plot(life_table$x, log10(life_table$lx/life_table$lx[1]))
points(x_vec,p_vec, col = "green")
points(x1,p_vec1, col = "green")
points(x2,p_vec2, col = "green")
points(x3,p_vec3, col = "green")
points(x4,p_vec4, col = "green")

age_vec = c("0-107","107-1079","1079-1216","1216-1771", "1771-1974")
mu_vec = c(result$par[1],result1$par[1],result2$par[1],result3$par[1],result4$par[1])
gamma_vec = c(result$par[2],result1$par[2],result2$par[2],result3$par[2],result4$par[2])

parameters = data.frame("Age subset" = age_vec, mu = mu_vec, gamma =gamma_vec, c = log10(mu_vec/(mu_vec + gamma_vec)))
