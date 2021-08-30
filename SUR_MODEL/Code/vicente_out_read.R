rm(list = ls())

output <- load("~/MAD_MODEL/VECTOR_MODEL/OUTPUT/param_VECTOR_MODEL1.RData")
out_mat <- matrix(0, nrow = 10, ncol = 2)
for(i in c(1:10)){
  out_mat[i,1] = lhs[[i]]$value
  out_mat[i,2] = lhs[[i]]$par
}

mat_sort <- out_mat[order(out_mat[,1], decreasing = F),]

