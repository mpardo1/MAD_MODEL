rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("coda")

# True param:
gam1 = 0.2
gam2 = 1.2
gam3 = 3
true1 = gam1
true2 = gam2
true3 = gam3
trueSD = 1
# Easy plots:
output <- load("~/Documentos/PHD/2021/SUR_Model/OUTPUT/chain_MH_2000eq_3param_5000_2021-09-21.RData")
# Ubuntu:
# output <- load("~/Documentos/PHD/2021/SUR_Model/RESULTS_ESTIMATION/MH/chain_MH_op_100eq_3param50000.RData")

# output2 <- load("~/Documents/PHD/2021/SUR_Model/PARAM_ESTIMATION/MH/Output/chain2_MH_op_3eq_3param1e+05.RData")

burnIn = 1000
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
### Summary: #######################

 par(mfrow = c(2,4))
 hist(chain[-(1:burnIn),1],nclass=30, main = expression(paste('Posterior of ', gamma[1])), xlab="True value = red line" )
 abline(v = mean(chain[-(1:burnIn),1]))
 abline(v = true1, col="red" )
 hist(chain[-(1:burnIn),2],nclass=30, main=expression(paste('Posterior of ', gamma[2])), xlab="True value = red line")
 abline(v = mean(chain[-(1:burnIn),2]))
 abline(v = true2, col="red" )
 hist(chain[-(1:burnIn),3],nclass=30, main=expression(paste('Posterior of ', gamma[3])), xlab="True value = red line")
 abline(v = mean(chain[-(1:burnIn),3]))#
 abline(v = true3, col="red" )
 hist(chain[-(1:burnIn),4],nclass=30, main="Posterior of sd", xlab="True value = red line")
 abline(v = mean(chain[-(1:burnIn),4]) )
 abline(v = trueSD, col="red" )

 plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = expression(paste('Chain values of ', gamma[1])), )
 abline(h = true1, col="red" )
 plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = expression(paste('Chain values of ', gamma[2])), )
 abline(h = true2, col="red" )
 plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = expression(paste('Chain values of ', gamma[3])), )
 abline(h = true3, col="red" )
 plot(chain[-(1:burnIn),4], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
 abline(h = trueSD, col="red" )

 
 # Plots from coda package:
 output <- load("~/Documentos/PHD/2021/SUR_Model/OUTPUT/chain_MH_op10000.RData")
 # Ubuntu:
 output <- load("~/Documentos/PHD/2021/SUR_Model/RESULTS_ESTIMATION/MH/chain_MH_op_100eq_3param50000.RData")
 chain_mc <- mcmc(chain)
 summary(chain_mc)
 plot(chain_mc)
 pairs(chain_mc)
 library(BayesianTools)
 correlationPlot(data.frame(chain_mc))
                 
  # Convergence diagnosis:
  rm(output)
  output2 <- load("~/Documents/PHD/2021/SUR_Model/OUTPUT/chain2_MH_op_100eq_3param10000.RData")
  # Ubuntu:
  output <- load("~/Documentos/PHD/2021/SUR_Model/OUTPUT/chain2_MH_op10000.RData")
  chain_mcmc2 <- mcmc(chain)
  plot(chain_mcmc2)
  combinedchains = mcmc.list(chain_mc, chain_mcmc2)
  plot(combinedchains)
  gelman.diag(combinedchains)
  gelman.plot(combinedchains)
 
  