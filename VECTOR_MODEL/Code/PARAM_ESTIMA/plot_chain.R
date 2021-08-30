rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("coda")

# True param:
true1 = 0.1
trueSD = 1
# Easy plots:
output <- load("~/Documents/PHD/2021/Mosquito_model/OUTPUT/chain2_MH_op10000.RData")
# output2 <- load("~/Documents/PHD/2021/SUR_Model/PARAM_ESTIMATION/MH/Output/chain2_MH_op_3eq_3param1e+05.RData")

burnIn = 500
acceptance = 1-mean(duplicated(chain[-(1:burnIn),]))
### Summary: #######################

par(mfrow = c(2,4))
hist(chain[-(1:burnIn),1],nclass=30, main="Posterior of a", xlab="True value = red line" )
abline(v = mean(chain[-(1:burnIn),1]))
abline(v = true1, col="red" )

hist(chain[-(1:burnIn),2],nclass=30, main="Posterior of sd", xlab="True value = red line")
abline(v = mean(chain[-(1:burnIn),2]) )
abline(v = trueSD, col="red" )

plot(chain[-(1:burnIn),1], type = "l", xlab="True value = red line" , main = "Chain values of a", )
abline(h = true1, col="red" )
plot(chain[-(1:burnIn),2], type = "l", xlab="True value = red line" , main = "Chain values of b", )
abline(h = true2, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of c", )
abline(h = true3, col="red" )
plot(chain[-(1:burnIn),3], type = "l", xlab="True value = red line" , main = "Chain values of sd", )
abline(h = trueSD, col="red" )


# Plots from coda package:
output <- load("~/Documents/PHD/2021/SUR_Model/PARAM_ESTIMATION/MH/Output/chain_MH_op_4eq_3param1e+05.RData")
chain_mc <- mcmc(chain)
summary(chain_mc)
plot(chain_mc)
pairs(chain_mc)
library(BayesianTools)
correlationPlot(data.frame(chain_mc))

# Convergence diagnosis:
rm(output)
output2 <- load("~/Documents/PHD/2021/SUR_Model/PARAM_ESTIMATION/MH/Output/chain2_MH_op_4eq_3param1e+05.RData")
chain_mcmc2 <- mcmc(chain)
plot(chain_mcmc2)
combinedchains = mcmc.list(chain_mc, chain_mcmc2)
plot(combinedchains)
gelman.diag(combinedchains)
gelman.plot(combinedchains)

