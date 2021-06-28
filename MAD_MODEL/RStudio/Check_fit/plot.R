xmin = 0 
xmax = 300
ggplot(int_sol) +
  geom_line(aes(X1,X2)) +
  geom_line(data = mat_fil, aes(x = X1, y = X2), color = "red",size = 0.6) +
  # xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 


ggplot(int_sol) +
  geom_line(aes(X1,X3)) +
  geom_line(data = mat_fil, aes(x = X1, y = X3), color = "red",size = 0.6) +
  # xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X4)) +
  geom_line(data = mat_fil, aes(x = X1, y = X4), color = "red",size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X5)) +
  geom_line(data = mat_fil, aes(x = X1, y = X5), color = "red",size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X6)) +
  geom_line(data = mat_fil, aes(x = X1, y = X6), color = "red",size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X32)) +
  geom_line(data = mat_fil, aes(x = X1, y = X32), color = "red",size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(mat_fil) +
  geom_line(aes(X1,X59), color = "red") +
  geom_line(data = int_sol, aes(x = X1, y = X59),size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(mat_fil) +
  geom_line(aes(X1,X80), color = "red") +
  geom_line(data = int_sol, aes(x = X1, y = X80),size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants")

ggplot(mat_fil) +
  geom_line(aes(X1,X16), color = "red") +
  geom_line(data = int_sol, aes(x = X1, y = X16),size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X16)) +
  geom_line(data = mat_fil, aes(x = X1, y = X16), color = "red",size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(mat_fil) +
  geom_line(aes(X1,X60), color = "red") +
  geom_line(data = int_sol, aes(x = X1, y = X60),size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X60)) +
  geom_line(data = mat_fil, aes(x = X1, y = X60), color = "red",size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X13)) +
  geom_line(data = mat_fil, aes(x = X1, y = X9), color = "red",size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X14)) +
  geom_line(data = mat_fil, aes(x = X1, y = X10), color = "red",size = 0.6) +
  # xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X15)) +
  geom_line(data = mat_fil, aes(x = X1, y = X11), color = "red",size = 0.6) +
  # xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X16)) +
  geom_line(data = mat_fil, aes(x = X1, y = X12), color = "red",size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

ggplot(int_sol) +
  geom_line(aes(X1,X17)) +
  geom_line(data = mat_fil, aes(x = X1, y = X13), color = "red",size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 


ggplot(int_sol) +
  geom_line(aes(X1,X25)) +
  geom_line(data = mat_fil, aes(x = X1, y = X21), color = "red",size = 0.6) +
  # xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 


ggplot(int_sol) +
  geom_line(aes(X1,X26)) +
  geom_line(data = mat_fil, aes(x = X1, y = X22), color = "red",size = 0.6) +
  #  xlim(xmin,xmax) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 


# Age = 30
ggplot(int_sol) +
  geom_line(aes(X1,X35)) +
  geom_line(data = mat_fil, aes(x = X1, y = X31), color = "red",size = 0.6) +
  #xlim(365,730) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

# Age = 231
ggplot(int_sol) +
  geom_line(aes(X1,X39)) +
  geom_line(data = mat_fil, aes(x = X1, y = X35), color = "red",size = 0.6) +
  # xlim(228,730) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 
# Age =  281
ggplot(int_sol) +
  geom_line(aes(X1,X40)) +
  geom_line(data = mat_fil, aes(x = X1, y = X36), color = "red",size = 0.6) +
  #  xlim(228,730) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 
# Age =  331
ggplot(int_sol) +
  geom_line(aes(X1,X41)) +
  geom_line(data = mat_fil, aes(x = X1, y = X37), color = "red",size = 0.6) +
  # xlim(330,730) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 
# Age = 1756
ggplot(int_sol) +
  geom_line(aes(X1,X60)) +
  geom_line(data = mat_fil, aes(x = X1, y = X56), color = "red",size = 0.6) +
  # xlim(1756,2000) + 
  xlab("Time (days)") + 
  ylab("Number of participants") 

