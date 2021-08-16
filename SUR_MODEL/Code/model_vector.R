rm(list = ls())
library(easypackages)
libraries("scales","gdata", "ggplot2", "numbers",
          "tidyverse","data.table","multiplex",
          "reshape","viridis","stats","ggpubr",
          "ggstatsplot","e1071","mlr3misc",
          "deSolve") 

# Registrations:
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")
head(down)

ggplot(down) + 
  geom_line(aes(x = time, y =down)) +
  ggtitle("Number of registrations in the app") +
  theme_bw()

min_t = 1
max_t = length(down$time)
for(i in c(min_t:(max_t-1))){
  dif = down$time[i+1] - down$time[i]
  if(dif > 1){
    print("gap")
    print(paste("i:",i))
  }
}

# Temperatures from Barcelona:
Path_temp = "~/MAD_MODEL/VECTOR_MODEL/data/bcn_weather_daily.Rds"
# Path_temp = paste(PC,Path_temp, sep="")

temp <-read_rds(Path_temp)
temp$date = as.Date(temp$date , "%Y-%m-%d")
temp <- temp %>%  group_by(date) %>% summarise(mean_temp = mean(valor))

ggplot(temp) +
  geom_line(aes(date, mean_temp))+
  ggtitle("Mean temperature Barcelona") + 
  xlab("Mean temperature")+
  theme_bw()

###############   ODE INTEGRATION   ##################
# require(deSolve)
# library.dynam.unload("deSolve", libpath=paste(.libPaths()[1], "//deSolve", sep=""))
# library.dynam("deSolve", package="deSolve", lib.loc=.libPaths()[1])
# OJOOOOO!!! Cuando cambias de PC borrar .o y .so.
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model.c")
dyn.load("model.so")

gam1 = 0.2
gam2 = 0.01
gam3 = 1.4
# We create a vector with the constant parameters.
parms = c(gam1,gam2,gam3)
# We set the initial conditions to cero.
# Y <- c(y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0, y5 = 0, y6 = 0, y7 = 0, y8 = 0, y9 = 0, y10 = 0,
#        y11 = 0.0, y12 = 0.0, y13 = 0.0, y14 = 0, y15 = 0, y16 = 0, y17 = 0, y18 = 0, y19 = 0, y20 = 0,
#        y21 = 0.0, y22 = 0.0, y23 = 0.0, y24 = 0, y25 = 0, y26 = 0, y27 = 0, y28 = 0, y29 = 0, y30 = 0,
#        y31 = 0.0, y32 = 0.0, y33 = 0.0, y34 = 0, y35 = 0, y36 = 0, y37 = 0, y38 = 0, y39 = 0, y40 = 0,
#        y41 = 0.0, y42 = 0.0, y43 = 0.0, y44 = 0, y45 = 0, y46 = 0, y47 = 0, y48 = 0, y49 = 0, y50 = 0,
#        y51 = 0.0, y52 = 0.0, y53 = 0.0, y54 = 0, y55 = 0, y56 = 0, y57 = 0, y58 = 0, y59 = 0, y60 = 0,
#        y61 = 0.0, y62 = 0.0, y63 = 0.0, y64 = 0, y65 = 0, y66 = 0, y67 = 0, y68 = 0, y69 = 0, y70 = 0,
#        y71 = 0.0, y72 = 0.0, y73 = 0.0, y74 = 0, y75 = 0, y76 = 0, y77 = 0, y78 = 0, y79 = 0, y80 = 0,
#        y81 = 0.0, y82 = 0.0, y83 = 0.0, y84 = 0, y85 = 0, y86 = 0, y87 = 0, y88 = 0, y89 = 0, y90 = 0,
#        y91 = 0.0, y92 = 0.0, y93 = 0.0, y94 = 0, y95 = 0, y96 = 0, y97 = 0, y98 = 0, y99 = 0, y100 = 0,
#        y111 = 0.0, y112 = 0.0, y113 = 0.0, y114 = 0, y115 = 0, y116 = 0, y117 = 0, y118 = 0, y119 = 0, y120 = 0,
#        y121 = 0.0, y122 = 0.0, y123 = 0.0, y124 = 0, y125 = 0, y126 = 0, y127 = 0, y128 = 0, y129 = 0, y130 = 0,
#        y131 = 0.0, y132 = 0.0, y133 = 0.0, y134 = 0, y135 = 0, y136 = 0, y137 = 0, y138 = 0, y139 = 0, y140 = 0,
#        y141 = 0.0, y142 = 0.0, y143 = 0.0, y144 = 0, y145 = 0, y146 = 0, y147 = 0, y148 = 0, y149 = 0, y150 = 0,
#        y151 = 0.0, y152 = 0.0, y153 = 0.0, y154 = 0, y155 = 0, y156 = 0, y157 = 0, y158 = 0, y159 = 0, y160 = 0,
#        y161 = 0.0, y162 = 0.0, y163 = 0.0, y164 = 0, y165 = 0, y166 = 0, y167 = 0, y168 = 0, y169 = 0, y170 = 0,
#        y171 = 0.0, y172 = 0.0, y173 = 0.0, y174 = 0, y175 = 0, y176 = 0, y177 = 0, y178 = 0, y179 = 0, y180 = 0,
#        y181 = 0.0, y182 = 0.0, y183 = 0.0, y184 = 0, y185 = 0, y186 = 0, y187 = 0, y188 = 0, y189 = 0, y190 = 0,
#        y191 = 0.0, y192 = 0.0, y193 = 0.0, y194 = 0, y195 = 0, y196 = 0, y197 = 0, y198 = 0, y199 = 0, y200 = 0,
#        y211 = 0.0, y212 = 0.0, y213 = 0.0, y214 = 0, y215 = 0, y216 = 0, y217 = 0, y218 = 0, y219 = 0, y220 = 0,
#        y221 = 0.0, y222 = 0.0, y223 = 0.0, y224 = 0, y225 = 0, y226 = 0, y227 = 0, y228 = 0, y229 = 0, y230 = 0,
#        y231 = 0.0, y232 = 0.0, y233 = 0.0, y234 = 0, y235 = 0, y236 = 0, y237 = 0, y238 = 0, y239 = 0, y240 = 0,
#        y241 = 0.0, y242 = 0.0, y243 = 0.0, y244 = 0, y245 = 0, y246 = 0, y247 = 0, y248 = 0, y249 = 0, y250 = 0,
#        y251 = 0.0, y252 = 0.0, y253 = 0.0, y254 = 0, y255 = 0, y256 = 0, y257 = 0, y258 = 0, y259 = 0, y260 = 0,
#        y261 = 0.0, y262 = 0.0, y263 = 0.0, y264 = 0, y265 = 0, y266 = 0, y267 = 0, y268 = 0, y269 = 0, y270 = 0,
#        y271 = 0.0, y272 = 0.0, y273 = 0.0, y274 = 0, y275 = 0, y276 = 0, y277 = 0, y278 = 0, y279 = 0, y280 = 0,
#        y281 = 0.0, y282 = 0.0, y283 = 0.0, y284 = 0, y285 = 0, y286 = 0, y287 = 0, y288 = 0, y289 = 0, y290 = 0,
#        y291 = 0.0, y292 = 0.0, y293 = 0.0, y294 = 0, y295 = 0, y296 = 0, y297 = 0, y298 = 0, y299 = 0, y300 = 0)
Y <- c(y1 = 0, y2 = 0, y3 = 0)

# List with the data frames of the forcings, sort as the c code.
forcs_mat <- list(data.matrix(down))

min_t <- min(down$time)
max_t <- max(down$time)
times <- seq(min_t,max_t, 1)
out <- ode(Y, times, func = "derivs",
           parms = parms, dllname = "model",
           initfunc = "initmod", nout = 1,
           outnames = "Sum", initforc = "forcc",
           forcings = down, 
           fcontrol = list(method = "constant")) 


ode <- data.frame(out)
ode$Sum <- NULL
head(ode)
df_plot <- reshape2::melt(ode, id.vars = c("time"))
unique(df_plot$variable)

df_plot_filt <- df_plot[df_plot$variable == 'y1' | df_plot$variable == 'y80',]

ggplot(df_plot_filt, aes(time, value))  +
  geom_line(aes( colour = variable)) +
  ylab("Counts") +
  ggtitle("Participation dynamics") +
  scale_color_manual(name = "",
                     labels = c("P1", "P80"),
                     values=c('#FF00F6','#FF2C00')) +
  theme_bw()

saveRDS(ode, file = "ode_pseudo.rds")

