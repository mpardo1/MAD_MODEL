rm(list = ls())
library(gdata) 
library(segmented)
library(e1071)
library(ggplot2)


# Read data from de csv and save it on life_table:
Path = "/home/marta/Documentos/PHD/2021/SUR_Model/John_Palmer/Life_table/participation_life_table.csv"
life_table = read.csv(Path) 

# Segmented does not allow referenced to the dataframe class
lx <- life_table$lx
x <- life_table$x
#First we create a linear model
LinearMod <- lm(lx~x)

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
num_breakp = 2
#seg_5 <- segmented(LinearMod, seg.Z = ~x, control=seg.control( fix.npsi=FALSE))
seg_5 <- segmented(LinearMod, seg.Z = ~x, npsi = num_breakp)
fit_data2 <- fitted(seg_5)
model <- data.frame(Age = x, Participants = fit_data2)
ggplot(model, aes(x = Age, y = Participants)) + geom_line() 
ggplot(model, aes(x = Age, y = Participants)) + geom_line() + geom_point(data = life_table, aes(x = x, y = lx), color = "red",legend.title = " nº breakpoints = 6")
ggplot(model, aes(x = Age, y = log10(Participants)))+ geom_line() + geom_point(data = life_table, aes(x = x, y = log10(lx)), color = "red")

# gamma_vec = c(1:num_breakp)
# # Create a function with the stationary state.
# n = num_breakp
# p_i <- function(i) {
#   gamma1 = gamma_vec[1]
#   gamma = gamma_vec[i]
#   p <- ((n-1)^(i-1) * gamma1)/((n)^(i-1)*gamma)
#   return(p)
# }
vec_psi = as.numeric(seg_5$psi[,2])[4:6]
# vec_psi[1] = vec_breakp[1]
# seg_mod <- segmented(LinearMod, seg.Z = ~x, seg.control(n.boot = 50,))
seg_mod <- segmented(LinearMod, seg.Z = ~x, psi = vec_psi)
fit_data2 <- fitted(seg_mod)
model <- data.frame(Age = x, Participants = fit_data2)
ggplot(model, aes(x = Age, y = Participants)) + geom_line() 
ggplot(model, aes(x = Age, y = Participants)) + geom_line() + geom_point(data = life_table, aes(x = x, y = lx), color = "red",legend.title = " nº breakpoints = 6")
ggplot(model, aes(x = Age, y = log10(Participants)))+ geom_line() + geom_point(data = life_table, aes(x = x, y = log10(lx)), color = "red")
