rm(list = ls())
# MAD model:
output <- load("~/MAD_MODEL/SUR_MODEL/Code/param_MAD_MODEL1.RData")
out_mat <- matrix(0, nrow = 10, ncol = 4)
for(i in c(1:10)){
  out_mat[i,1] = lhs[[i]]$value
  out_mat[i,2:4] = lhs[[i]]$par
}

mat_sort <- out_mat[order(out_mat[,1], decreasing = F),]

# Vector model:
output <- load("~/Documents/PHD/2021/Mosquito_model/OUTPUT/NM/param_VECTOR_MODEL1719.RData")
out_mat <- matrix(0, nrow = 1000, ncol = 2)
for(i in c(1:10)){
  out_mat[i,1] = lhs[[i]]$value
  out_mat[i,2] = lhs[[i]]$par
}

mat_sort <- out_mat[order(out_mat[,1], decreasing = F),]

