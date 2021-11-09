rm(list = ls())
library(gdata) 
library(segmented)
library(e1071)
library(ggplot2)


# Read data from de csv and save it on life_table:
Path = "~/MAD_MODEL/SUR_MODEL/data/participation_life_table.csv"
life_table = read.csv(Path) 

# Segmented does not allow referenced to the dataframe class
lx <- life_table$lx
x <- life_table$x
#First we create a linear model
LinearMod <- lm(lx~x)

#Create vectors to save the squared errors and the number of breakpoints.
rSq_ve = c(1:8)
nod_vec = c(1:8)

# Create a loop over the fixed number of breakpoints for the segmented algorithm.
for(i in 1:8){
  seg <- segmented(LinearMod, seg.Z = ~x, npsi = i)
  nod_vec[i] = i
  rSq_ve[i] = summary(seg)$adj.r.squared
}

df_break <- data.frame(nod = nod_vec, r_sq = rSq_ve)
plot_J <- ggplot(df_break) + geom_line(aes(nod,r_sq)) 
plot_J + theme_bw() + xlab("Number of breakpoints") + ylab("R squared")
#Plot the squared error against the number of breakpoints.
plot(nod_vec, rSq_ve, xlab = "Number of breakpoints", ylab = "R squared", xlim = c(0,10))

# Do segmented for 5 breakpoints since it appears to be the most efficient number of breakpoints.
num_breakp = 2
#seg_5 <- segmented(LinearMod, seg.Z = ~x, control=seg.control( fix.npsi=FALSE))
seg_5 <- segmented(LinearMod, seg.Z = ~x, npsi = num_breakp)
fit_data2 <- fitted(seg_5)
model <- data.frame(Age = x, estimation = fit_data2, lx = life_table$lx)
df_plot <- reshape2::melt(model, id.vars = c("Age"))

ggplot(life_table, aes(x = x, y = lx)) +
  geom_line() + theme_bw() +
  xlab("Age of the participants")

ggplot(df_plot,aes(Age, value)) +
  geom_line(aes( colour = variable)) + theme_bw() +
  scale_color_manual(name = "",
                     values=c('#3c33ff','#ff3933'))+
  xlab("Age of the participants")

model <- data.frame(Age = x, estimation = log10(fit_data2), lx = log10(life_table$lx))
df_plot <- reshape2::melt(model, id.vars = c("Age"))
ggplot(df_plot,aes(Age, value)) +
  geom_line(aes( colour = variable)) + theme_bw() +
  scale_color_manual(name = "",
                     values=c('#3c33ff','#ff3933'))+
  xlab("Age of the participants")

vec_psi = as.numeric(seg_5$psi[,2])[4:6]
# vec_psi[1] = vec_breakp[1]
# seg_mod <- segmented(LinearMod, seg.Z = ~x, seg.control(n.boot = 50,))
seg_mod <- segmented(LinearMod, seg.Z = ~x, psi = vec_psi)
fit_data2 <- fitted(seg_mod)
model <- data.frame(Age = x, Participants = fit_data2)
ggplot(model, aes(x = Age, y = Participants)) + geom_line() 
ggplot(model, aes(x = Age, y = Participants)) + geom_line() + geom_point(data = life_table, aes(x = x, y = lx), color = "red",legend.title = " nº breakpoints = 6")
ggplot(model, aes(x = Age, y = log10(Participants)))+ geom_line() + geom_point(data = life_table, aes(x = x, y = log10(lx)), color = "red")

