rm(list = ls())
library(gdata) 
library(segmented)
library(e1071)
library(ggplot2)


# Read data from de csv and save it on life_table:
Path = "~/MAD_MODEL/SUR_MODEL/data/participation_life_table.csv"
life_table = read.csv(Path) 

# Quadratic fit:
life_table_filt <- life_table[2:177,]
ggplot(life_table_filt) +
  geom_point(aes(x, nqx))+
  theme_bw()

y <- life_table_filt$nqx
x <- life_table_filt$x
x2 <- x*x
model <- glm(y ~ x + x2)

#plot x vs. y
plot(x, y)

#define x-values to use for regression line
l_out <- 10000
x=seq(from=0,to=1200,length.out=l_out)

#use the model to predict the y-values based on the x-values
y=predict(model,list(days=x, x2=x^2))

#add the fitted regression line to the plot (lwd specifies the width of the line)
matlines(x,y, lwd=2)

df_log <- data.frame(x,y)
df_log <- df_log[-1,]
ggplot(life_table_filt) + 
  geom_point(aes(x,nqx), size = 1) +
  geom_line(data=df_log, aes(x,y), color = "#863390") + 
  theme_bw() + xlab("Days since registration") +
  ylab("Weekly mortality (probability)")

# Quadratic fit lx:
life_table_filt <- life_table[2:177,]
ggplot(life_table_filt) +
  geom_point(aes(x, lx))+
  theme_bw()

y <- life_table_filt$lx
x <- life_table_filt$x
model <- glm(y ~ log(x))

#plot x vs. y
plot(x, y)

#define x-values to use for regression line
l_out <- 10000
x=seq(from=0,to=1200,length.out=l_out)

#use the model to predict the y-values based on the x-values
y=predict(model,list(days=x))

#add the fitted regression line to the plot (lwd specifies the width of the line)
matlines(x,y, lwd=2)

df_log <- data.frame(x,y)
df_log <- df_log[-1,]
ggplot(life_table_filt) + 
  geom_point(aes(x,lx), size = 1) +
  # geom_line(data=df_log, aes(x,y), color = "#863390") + 
  theme_bw() + xlab("Days since registration") +
  ylab("Number of participants")

# Geometric mean:
max_x <-  nrow(life_table)
x_geo <- c()
init <- c()
x_geo[1] <- 0
i <- 2
tol <-  0
init[1] <- 0
while(tol < max_x){
  x_geo[i] <- x_geo[i-1] + 2^(i-2)
  tol <- x_geo[i] + 2^(i-1)
  # init[i] <- which(life_table$x > x_geo[i] & life_table$x < tol )
  i <-  i + 1
}
# Because it starts in 0:
x_geo <- x_geo +  1 

# Function  to compute the geometric mean:
geo_mean <- function(x){exp(mean(log(x)))} 
mean_g <- matrix(0,length(x_geo), 4 )
for(i in c(1:(length(x_geo)-1))){
  init <- x_geo[i] 
  end <- x_geo[i+1]
  mean_g[i,] <- c(life_table$x[init], life_table$x[end],
                geo_mean(life_table$x[init:end]),
                geo_mean(life_table$nqx[init:end]))
}

df_mg <- as.data.frame(mean_g)
df_mg <-  df_mg[-nrow(df_mg),]
colnames(df_mg) <-  c("left", "right", "g_mean_int", "g_mean_prob")
ggplot(df_mg) + 
  geom_point(aes(g_mean_int,g_mean_prob), size = 0.8)+
  theme_bw() + xlab("Days since registration") +
  ylab("Mortality probability")

y <- df_mg$g_mean_prob
x <- df_mg$g_mean_int
model <- lm(y ~ log(x))

#plot x vs. y
plot(x, y)

#define x-values to use for regression line
l_out <- 10000
x=seq(from=0,to=1500,length.out=l_out)

#use the model to predict the y-values based on the x-values
y=predict(model,newdata=list(x=seq(from=0,to=1500,length.out=l_out)),
          interval="confidence")

#add the fitted regression line to the plot (lwd specifies the width of the line)
matlines(x,y, lwd=2)

df_log <- data.frame(x,y)
df_log <- df_log[-1,]
ggplot(df_mg) + 
  geom_point(aes(g_mean_int,g_mean_prob), size = 1) +
  geom_line(data=df_log, aes(x,fit), color = "#863390") + 
  theme_bw() + xlab("Days since registration") +
  ylab("Mortality probability")
# nqx: prob of dying within the next n=7 days.
nqx_comp <- (life_table$lx - life_table$lxn)/life_table$lx
df_qx <- data.frame(age = life_table$x, nqx = life_table$nqx, nqx_comp)

ggplot(df_qx) +
  geom_line(aes(age, nqx))+
  theme_bw()

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
  # xlab("Age of the participants") +
  ylab("People \"surviving\"")

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

