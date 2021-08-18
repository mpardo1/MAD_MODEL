rm(list = ls())

output <- load("~/Documents/chain_MH_1e+05.RData")
out_mat <- matrix(0, nrow = 1000, ncol = 4)
for(i in c(1:300)){
  out_mat[i,1] = lhs[[i]]$value
  out_mat[i,2:4] = lhs[[i]]$par
}

mat_sort <- out_mat[order(out_mat[,1], decreasing = F),]
