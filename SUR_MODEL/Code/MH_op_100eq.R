rm(list = ls())
start_time <- Sys.time()

library("parallel")
library("tidyverse")
library("deSolve")
library("coda")
# Create Pseudo data:
# Create Pseudo Data:
Path = "~/MAD_MODEL/SUR_MODEL/Code/"
# Path = paste(PC,Path, sep="")

setwd(Path)
system("R CMD SHLIB model.c")
dyn.load("model.so")

gam1 = 0.2
gam2 = 1.2
gam3 = 3
# We create a vector with the constant parameters.
parms = c(gam1,gam2,gam3)
# We set the initial conditions to zero.
Y <- c(y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0, y5 = 0, y6 = 0, y7 = 0, y8 = 0, y9 = 0, y10 = 0,
       y11 = 0.0, y12 = 0.0, y13 = 0.0, y14 = 0, y15 = 0, y16 = 0, y17 = 0, y18 = 0, y19 = 0, y20 = 0,
       y21 = 0.0, y22 = 0.0, y23 = 0.0, y24 = 0, y25 = 0, y26 = 0, y27 = 0, y28 = 0, y29 = 0, y30 = 0,
       y31 = 0.0, y32 = 0.0, y33 = 0.0, y34 = 0, y35 = 0, y36 = 0, y37 = 0, y38 = 0, y39 = 0, y40 = 0,
       y41 = 0.0, y42 = 0.0, y43 = 0.0, y44 = 0, y45 = 0, y46 = 0, y47 = 0, y48 = 0, y49 = 0, y50 = 0,
       y51 = 0.0, y52 = 0.0, y53 = 0.0, y54 = 0, y55 = 0, y56 = 0, y57 = 0, y58 = 0, y59 = 0, y60 = 0,
       y61 = 0.0, y62 = 0.0, y63 = 0.0, y64 = 0, y65 = 0, y66 = 0, y67 = 0, y68 = 0, y69 = 0, y70 = 0,
       y71 = 0.0, y72 = 0.0, y73 = 0.0, y74 = 0, y75 = 0, y76 = 0, y77 = 0, y78 = 0, y79 = 0, y80 = 0,
       y81 = 0.0, y82 = 0.0, y83 = 0.0, y84 = 0, y85 = 0, y86 = 0, y87 = 0, y88 = 0, y89 = 0, y90 = 0,
       y91 = 0.0, y92 = 0.0, y93 = 0.0, y94 = 0, y95 = 0, y96 = 0, y97 = 0, y98 = 0, y99 = 0, y100 = 0)

Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")

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

saveRDS(ode, file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo.rds")
# Función de c que corre la ODE -------------------------------------------

# Esta función está adaptada para su uso con el paquete deSolve, el estándar en
# R para las ODE
true1 = gam1
true2 = gam2
true3 = gam3
trueSD = 1

# Likelihood;
likelihood <- function(x) # forzamientos para el solver de la ode
{ 
  if(x[1] < 0 | x[2] < 0 | x[3] < 0 ){
    res = -86829146000
  }else{
    pars <- c(gam1 = x[1], # death rate group 1
              gam2 = x[2],
              gam3 = x[3]) # death rate group 2
    
    sd <- x[4]
    
    population <- c(y1 = 0.0, y2 = 0.0, y3 = 0.0, y4 = 0, y5 = 0, y6 = 0, y7 = 0, y8 = 0, y9 = 0, y10 = 0,
                    y11 = 0.0, y12 = 0.0, y13 = 0.0, y14 = 0, y15 = 0, y16 = 0, y17 = 0, y18 = 0, y19 = 0, y20 = 0,
                    y21 = 0.0, y22 = 0.0, y23 = 0.0, y24 = 0, y25 = 0, y26 = 0, y27 = 0, y28 = 0, y29 = 0, y30 = 0,
                    y31 = 0.0, y32 = 0.0, y33 = 0.0, y34 = 0, y35 = 0, y36 = 0, y37 = 0, y38 = 0, y39 = 0, y40 = 0,
                    y41 = 0.0, y42 = 0.0, y43 = 0.0, y44 = 0, y45 = 0, y46 = 0, y47 = 0, y48 = 0, y49 = 0, y50 = 0,
                    y51 = 0.0, y52 = 0.0, y53 = 0.0, y54 = 0, y55 = 0, y56 = 0, y57 = 0, y58 = 0, y59 = 0, y60 = 0,
                    y61 = 0.0, y62 = 0.0, y63 = 0.0, y64 = 0, y65 = 0, y66 = 0, y67 = 0, y68 = 0, y69 = 0, y70 = 0,
                    y71 = 0.0, y72 = 0.0, y73 = 0.0, y74 = 0, y75 = 0, y76 = 0, y77 = 0, y78 = 0, y79 = 0, y80 = 0,
                    y81 = 0.0, y82 = 0.0, y83 = 0.0, y84 = 0, y85 = 0, y86 = 0, y87 = 0, y88 = 0, y89 = 0, y90 = 0,
                    y91 = 0.0, y92 = 0.0, y93 = 0.0, y94 = 0, y95 = 0, y96 = 0, y97 = 0, y98 = 0, y99 = 0, y100 = 0)
    #Vector inicial para ODE
    
    forcs_mat <- list(data.matrix(forcings))
    
    z <- ode(y = population,
             times = 0:nrow(y), func = "derivs", method = "ode45",
             dllname = "model", initfunc = "initmod", nout = 0, 
             parms = pars, initforc = "forcc", forcings = forcs_mat, 
             fcontrol = list(method = "constant")) #Aquí corre el ODE
    
    colnames(z)[2:101] <- c("P1", "P2", "P3", "P4", "P5","P6", "P7", "P8", "P9", "P10",
                            "P11", "P12", "P13", "P14", "P15","P16", "P17", "P18", "P19", "P20",
                            "P21", "P22", "P23", "P24", "P25","P26", "P27", "P28", "P29", "P30",
                            "P31", "P32", "P33", "P34", "P35","P36", "P37", "P38", "P39", "P40",
                            "P41", "P42", "P43", "P44", "P45","P46", "P47", "P48", "P49", "P50",
                            "P51", "P52", "P53", "P54", "P55","P56", "P57", "P58", "P59", "P60",
                            "P61", "P62", "P63", "P64", "P65","P66", "P67", "P68", "P69", "P70",
                            "P71", "P72", "P73", "P74", "P75","P76", "P77", "P78", "P79", "P80",
                            "P81", "P82", "P83", "P84", "P85","P86", "P87", "P88", "P89", "P90",
                            "P91", "P92", "P93", "P94", "P95","P96", "P97", "P98", "P99", "P100"
    )
    
    z <- as.data.frame(z)
    z <- z[-1, ]
    
    P1 <- y$X1;  P2 <- y$X2;  P3 <- y$X3;  P4 <-  y$X4; P5 <- y$X5;  P6 <- y$X6  ;P7 <- y$X7;  P8 <- y$X8;   P9 <- y$X9;P10 <- y$X10
    P11 <- y$X11;P12 <- y$X12;P13 <- y$X13;P14 <- y$X14;P15 <- y$X15;P16 <- y$X16;P17 <- y$X17;P18 <- y$X18; P19 <- y$X19;P20 <- y$X20
    P21 <- y$X21;P22 <- y$X22;P23 <- y$X23;P24 <- y$X24;P25 <- y$X25;P26 <- y$X26;P27 <- y$X27;P28 <- y$X28; P29 <- y$X29;P30 <- y$X30
    P31 <- y$X31;P32 <- y$X32;P33 <- y$X33;P34 <- y$X34;P35 <- y$X35;P36 <- y$X36;P37 <- y$X37;P38 <- y$X38; P39 <- y$X39;P40 <- y$X40
    P41 <- y$X41;P42 <- y$X42;P43 <- y$X43;P44 <- y$X44;P45 <- y$X45;P46 <- y$X46;P47 <- y$X47;P48 <- y$X48; P49 <- y$X49;P50 <- y$X50
    P51 <- y$X51;P52 <- y$X52;P53 <- y$X53;P54 <- y$X54;P55 <- y$X55;P56 <- y$X56;P57 <- y$X57;P58 <- y$X58; P59 <- y$X59;P60 <- y$X60
    P61 <- y$X61;P62 <- y$X62;P63 <- y$X63;P64 <- y$X64;P65 <- y$X65;P66 <- y$X66;P67 <- y$X67;P68 <- y$X68; P69 <- y$X69;P70 <- y$X70
    P71 <- y$X71;P72 <- y$X72;P73 <- y$X73;P74 <- y$X74;P75 <- y$X75;P76 <- y$X76;P77 <- y$X77;P78 <- y$X78; P79 <- y$X79;P80 <- y$X80
    P81 <- y$X81;P82 <- y$X82;P83 <- y$X83;P84 <- y$X84;P85 <- y$X85;P86 <- y$X86;P87 <- y$X87;P88 <- y$X88; P89 <- y$X89;P90 <- y$X90
    P91 <- y$X91;P92 <- y$X92;P93 <- y$X93;P94 <- y$X94;P95 <- y$X95;P96 <- y$X96;P97 <- y$X97;P98 <- y$X98; P99 <- y$X99;P100 <- y$X100
    
    
    res <- #cálculo de la loglikelihood en función de las desviaciones estándar
      sum(dnorm(P1, mean = z$P1, sd = sd, log = T)) + sum(dnorm(P6, mean = z$P6, sd = sd, log = T)) +
      sum(dnorm(P2, mean = z$P2, sd = sd, log = T)) + sum(dnorm(P7, mean = z$P7, sd = sd, log = T)) +
      sum(dnorm(P3, mean = z$P3, sd = sd, log = T)) + sum(dnorm(P8, mean = z$P8, sd = sd, log = T)) +
      sum(dnorm(P4, mean = z$P4, sd = sd, log = T)) + sum(dnorm(P9, mean = z$P9, sd = sd, log = T)) +
      sum(dnorm(P5, mean = z$P5, sd = sd, log = T)) + sum(dnorm(P10, mean = z$P10, sd = sd, log = T)) +
      sum(dnorm(P11, mean = z$P11, sd = sd, log = T)) + sum(dnorm(P16, mean = z$P16, sd = sd, log = T)) +
      sum(dnorm(P12, mean = z$P12, sd = sd, log = T)) + sum(dnorm(P17, mean = z$P17, sd = sd, log = T)) +
      sum(dnorm(P13, mean = z$P13, sd = sd, log = T)) + sum(dnorm(P18, mean = z$P18, sd = sd, log = T)) +
      sum(dnorm(P14, mean = z$P14, sd = sd, log = T)) + sum(dnorm(P19, mean = z$P19, sd = sd, log = T)) +
      sum(dnorm(P15, mean = z$P15, sd = sd, log = T)) + sum(dnorm(P20, mean = z$P20, sd = sd, log = T)) +
      sum(dnorm(P21, mean = z$P21, sd = sd, log = T)) + sum(dnorm(P26, mean = z$P26, sd = sd, log = T)) +
      sum(dnorm(P22, mean = z$P22, sd = sd, log = T)) + sum(dnorm(P27, mean = z$P27, sd = sd, log = T)) +
      sum(dnorm(P23, mean = z$P23, sd = sd, log = T)) + sum(dnorm(P28, mean = z$P28, sd = sd, log = T)) +
      sum(dnorm(P24, mean = z$P24, sd = sd, log = T)) + sum(dnorm(P29, mean = z$P29, sd = sd, log = T)) +
      sum(dnorm(P25, mean = z$P25, sd = sd, log = T)) + sum(dnorm(P30, mean = z$P30, sd = sd, log = T)) +
      sum(dnorm(P31, mean = z$P31, sd = sd, log = T)) + sum(dnorm(P36, mean = z$P36, sd = sd, log = T)) +
      sum(dnorm(P32, mean = z$P32, sd = sd, log = T)) + sum(dnorm(P37, mean = z$P37, sd = sd, log = T)) +
      sum(dnorm(P33, mean = z$P33, sd = sd, log = T)) + sum(dnorm(P38, mean = z$P38, sd = sd, log = T)) +
      sum(dnorm(P34, mean = z$P34, sd = sd, log = T)) + sum(dnorm(P39, mean = z$P39, sd = sd, log = T)) +
      sum(dnorm(P35, mean = z$P35, sd = sd, log = T)) + sum(dnorm(P40, mean = z$P40, sd = sd, log = T)) +
      sum(dnorm(P41, mean = z$P41, sd = sd, log = T)) + sum(dnorm(P46, mean = z$P46, sd = sd, log = T)) +
      sum(dnorm(P42, mean = z$P42, sd = sd, log = T)) + sum(dnorm(P47, mean = z$P47, sd = sd, log = T)) +
      sum(dnorm(P43, mean = z$P43, sd = sd, log = T)) + sum(dnorm(P48, mean = z$P48, sd = sd, log = T)) +
      sum(dnorm(P44, mean = z$P44, sd = sd, log = T)) + sum(dnorm(P49, mean = z$P49, sd = sd, log = T)) +
      sum(dnorm(P45, mean = z$P45, sd = sd, log = T)) + sum(dnorm(P50, mean = z$P50, sd = sd, log = T)) +
      sum(dnorm(P51, mean = z$P51, sd = sd, log = T)) + sum(dnorm(P56, mean = z$P56, sd = sd, log = T)) +
      sum(dnorm(P52, mean = z$P52, sd = sd, log = T)) + sum(dnorm(P57, mean = z$P57, sd = sd, log = T)) +
      sum(dnorm(P53, mean = z$P53, sd = sd, log = T)) + sum(dnorm(P58, mean = z$P58, sd = sd, log = T)) +
      sum(dnorm(P54, mean = z$P54, sd = sd, log = T)) + sum(dnorm(P59, mean = z$P59, sd = sd, log = T)) +
      sum(dnorm(P55, mean = z$P55, sd = sd, log = T)) + sum(dnorm(P60, mean = z$P60, sd = sd, log = T)) +
      sum(dnorm(P61, mean = z$P61, sd = sd, log = T)) + sum(dnorm(P66, mean = z$P66, sd = sd, log = T)) +
      sum(dnorm(P62, mean = z$P62, sd = sd, log = T)) + sum(dnorm(P67, mean = z$P67, sd = sd, log = T)) +
      sum(dnorm(P63, mean = z$P63, sd = sd, log = T)) + sum(dnorm(P68, mean = z$P68, sd = sd, log = T)) +
      sum(dnorm(P64, mean = z$P64, sd = sd, log = T)) + sum(dnorm(P69, mean = z$P69, sd = sd, log = T)) +
      sum(dnorm(P65, mean = z$P65, sd = sd, log = T)) + sum(dnorm(P70, mean = z$P70, sd = sd, log = T)) +
      sum(dnorm(P71, mean = z$P71, sd = sd, log = T)) + sum(dnorm(P76, mean = z$P76, sd = sd, log = T)) +
      sum(dnorm(P72, mean = z$P72, sd = sd, log = T)) + sum(dnorm(P77, mean = z$P77, sd = sd, log = T)) +
      sum(dnorm(P73, mean = z$P73, sd = sd, log = T)) + sum(dnorm(P78, mean = z$P78, sd = sd, log = T)) +
      sum(dnorm(P74, mean = z$P74, sd = sd, log = T)) + sum(dnorm(P79, mean = z$P79, sd = sd, log = T)) +
      sum(dnorm(P75, mean = z$P75, sd = sd, log = T)) + sum(dnorm(P80, mean = z$P80, sd = sd, log = T)) +
      sum(dnorm(P81, mean = z$P81, sd = sd, log = T)) + sum(dnorm(P86, mean = z$P86, sd = sd, log = T)) +
      sum(dnorm(P82, mean = z$P82, sd = sd, log = T)) + sum(dnorm(P87, mean = z$P87, sd = sd, log = T)) +
      sum(dnorm(P83, mean = z$P83, sd = sd, log = T)) + sum(dnorm(P88, mean = z$P88, sd = sd, log = T)) +
      sum(dnorm(P84, mean = z$P84, sd = sd, log = T)) + sum(dnorm(P89, mean = z$P89, sd = sd, log = T)) +
      sum(dnorm(P85, mean = z$P85, sd = sd, log = T)) + sum(dnorm(P90, mean = z$P90, sd = sd, log = T)) +
      sum(dnorm(P91, mean = z$P91, sd = sd, log = T)) + sum(dnorm(P96, mean = z$P96, sd = sd, log = T)) +
      sum(dnorm(P92, mean = z$P92, sd = sd, log = T)) + sum(dnorm(P97, mean = z$P97, sd = sd, log = T)) +
      sum(dnorm(P93, mean = z$P93, sd = sd, log = T)) + sum(dnorm(P98, mean = z$P98, sd = sd, log = T)) +
      sum(dnorm(P94, mean = z$P94, sd = sd, log = T)) + sum(dnorm(P99, mean = z$P99, sd = sd, log = T)) +
      sum(dnorm(P95, mean = z$P95, sd = sd, log = T)) + sum(dnorm(P100, mean = z$P100, sd = sd, log = T)) 
  }
  return(res)
}


# Registrations:
Path = "~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data"
down <- data.frame(t(read.table(Path, header=FALSE)))
colnames(down) <- c("time", "down")
head(down)
forcs_mat <- data.matrix(down)

# Pseudo Data to check the oprimization method.
ob_data <- readRDS(file = "~/MAD_MODEL/SUR_MODEL/Code/ode_pseudo.rds")
colnames(ob_data) <- c("time", "X1", "X2", "X3", "X4", "X5" , "X6", "X7", "X8", "X9", "X10",
                       "X11", "X12", "X13", "X14", "X15" , "X16", "X17", "X18", "X19", "X20",
                       "X21", "X22", "X23", "X24", "X25" , "X26", "X27", "X28", "X29", "X30",
                       "X31", "X32", "X33", "X34", "X35" , "X36", "X37", "X38", "X39", "X40",
                       "X41", "X42", "X43", "X44", "X45" , "X46", "X47", "X48", "X49", "X50",
                       "X51", "X52", "X53", "X54", "X55" , "X56", "X57", "X58", "X59", "X60",
                       "X61", "X62", "X63", "X64", "X65" , "X66", "X67", "X68", "X69", "X70",
                       "X71", "X72", "X73", "X74", "X75" , "X76", "X77", "X78", "X79", "X80",
                       "X81", "X82", "X83", "X84", "X85" , "X86", "X87", "X88", "X89", "X90",
                       "X91", "X92", "X93", "X94", "X95" , "X96", "X97", "X98", "X99", "X100")

l <- nrow(ob_data)
ob_data[,2:101] <- ob_data[,2:101] + matrix(rnorm(l*100,0,trueSD), ncol = 100, nrow = l)


y <- ob_data
forcings <- down
# head(ob_data)

# Prior distribution
prior = function(param){
  a = param[1]
  b = param[2]
  c = param[3]
  sd = param[4]
  aprior = dnorm(a, sd=1,  log = T)
  bprior = dnorm(b, sd=1,  log = T)
  cprior = dnorm(c, sd=1,  log = T)
  sdprior = dnorm(sd, sd=1,  log = T)
  return(aprior+bprior+cprior+sdprior)
}

######## Metropolis algorithm ################

proposalfunction = function(param){
  vec <- param + c(rnorm(3, mean = c(0,0,0), sd= c(0.1,0.1,0.1))
  ,abs(rnorm(1,mean = 0 ,sd = 0.3)))
  # vec <- c(rnorm(3, mean = param[1:3], sd= c(0.1,0.1,0.1))
  #          ,abs(rnorm(1,mean = param[4] ,sd = 0.1)))
  
  return(vec)
}

run_metropolis_MCMC = function(startvalue, iterations){
  chain = array(dim = c(iterations+1,4))
  chain[1,] = startvalue
  for (i in 1:iterations){
    proposal = proposalfunction(chain[i,])
    print("Iteration:")
    print(i)
    probab = exp(likelihood(proposal)+ prior(proposal) - likelihood(chain[i,])- prior(chain[i,]))
    if (runif(1) < probab){
      chain[i+1,] = proposal
    }else{
      chain[i+1,] = chain[i,] } 
  }
  return(chain)
}


startvalue = c(0.1,1,2.5,0.5)
iterations = 10000
chain = run_metropolis_MCMC(startvalue, iterations)

filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/chain_MH_op_100eq_3param",iterations,".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain, file = filename)
# burnIn = 5000
# acceptance1 = 1-mean(duplicated(chain[-(1:burnIn),]))
# 
# chain <- mcmc(chain)
# summary(chain)
# plot(chain)

# library(BayesianTools)
# correlationPlot(data.frame(chain))
# #
# # Convergence diagnosis:
#
# print("Optimization finish")
chain2 = run_metropolis_MCMC(startvalue, iterations)
# burnIn = 5000
# acceptance2 = 1-mean(duplicated(chain2[-(1:burnIn),]))

# chain2 <- mcmc(chain2)
filename <- paste0("~/MAD_MODEL/SUR_MODEL/Code/chain2_MH_op_100eq_3param",iterations,".RData") #Salva cada ronda de optimizaciones, por si acaso
save(chain2, file = filename)

# 
# combinedchains = mcmc.list(chain, chain2)
# plot(combinedchains)
# gelman.diag(combinedchains)
# gelman.plot(combinedchains)
end_time <- Sys.time()
diff_time <- end_time - start_time
print("Execution time:")
print(diff_time)