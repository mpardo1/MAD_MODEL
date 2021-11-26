rm(list = ls())
library("parallel")
library("tidyverse")
library("deSolve")
library("ggpubr")
library("ggsci")
library("gridExtra")
library("grid")
library("kableExtra")

# Read output from param_estimation.R
output <- load("/home/marta/Documentos/PHD/2021/SUR_Model/OUTPUT/param_MAD_MODEL_1core_900it1.RData")
output <- load("~/Documents/PHD/2021/SUR_Model/OUTPUT/param_MAD_MODEL_1core_900it1.RData")

# Number of seeds:
N = 300
out_mat <- matrix(0, nrow = N, ncol = 4)
for(i in c(1:N)){
  out_mat[i,1] = lhs[[i]]$value
  out_mat[i,2:4] = lhs[[i]]$par
}

mat_sort <- out_mat[order(out_mat[,1], decreasing = T),]

# Load C code:
Path = "~/MAD_MODEL/SUR_MODEL/Code/PARAM_ESTIMATION/MH/"
setwd(Path)
system("R CMD SHLIB model_2000eq.c")
dyn.load("model_2000eq.so")

# --------------------- NEW DATA ------------------------
# Read Observed data:
ob_data <- read.table("~/MAD_MODEL/SUR_MODEL/data/Observed_2021-11-23_2691.dat", header=FALSE, sep= "\t")
ob_data_t <- t(ob_data)
downloads <- ob_data_t[,1:2]
end_time <- max(ob_data_t[,1])
age_max <- max(ob_data_t[1,]) -1

# ------------------ OLD DATA ----------------------------
ob_data <- read.table("~/MAD_MODEL/SUR_MODEL/data/Observed_data_2300.data", header=FALSE, sep= " ")
ob_data_t <- t(ob_data)
downloads <- t(read.table("~/MAD_MODEL/SUR_MODEL/data/Downloads_2378.data", header=FALSE, sep= " "))
end_time <- max(ob_data[1,])

# ------------------INTEGRATION --------------------------
# Parameters:
for(i in c(1:10)){
    it <- i
    gam1 = mat_sort[it,2]
    gam2 = mat_sort[it,3]
    gam3 = mat_sort[it,4]
    
    pars <- c(gam1,gam2,gam3)
    population <- matrix(0,nrow = 1,ncol=2000)
    z <- ode(y = population,
             times = 0:end_time, func = "derivs", method = "ode45",
             dllname = "model_2000eq", initfunc = "initmod", nout = 0, 
             parms = pars, initforc = "forcc", forcings = downloads, 
             fcontrol = list(method = "constant"))
    
    z <- as.data.frame(z)
    
    # ----------------- 1 DAY OLD --------------------------#
    df_age1 <- data.frame( time = z$time, ob = ob_data_t[,2], sim = z$`1`) 
    df_plot <- reshape2::melt(df_age1, id.vars = c("time"))
    
    # Plot number of participants at each time.
    plot_1  <- ggplot(df_plot,aes(time, value)) + 
      geom_line(aes( colour = variable))  +
      ylab("Number of participants") + 
      theme_bw() + scale_color_nejm()
    
    df_age30 <- data.frame( time = z$time, ob = ob_data_t[,31], sim = z$`30`) 
    df_plot <- reshape2::melt(df_age30, id.vars = c("time"))
    
    # Plot number of participants at each time.
    plot_31  <-  ggplot(df_plot,aes(time, value)) + 
      geom_line(aes( colour = variable))  +
      ylab("Number of participants") +
      scale_color_nejm() + 
      theme_bw()
    
    df_age300 <- data.frame( time = z$time, ob = ob_data_t[,301], sim = z$`300`) 
    df_plot <- reshape2::melt(df_age300, id.vars = c("time"))
    
    # Plot number of participants at each time.
    plot_300  <- ggplot(df_plot,aes(time, value)) + 
      geom_line(aes( colour = variable))  +
      ylab("Number of participants") +
      scale_color_nejm() + 
      theme_bw()
    
    df_age600 <- data.frame( time = z$time, ob = ob_data_t[,601], sim = z$`600`) 
    df_plot <- reshape2::melt(df_age600, id.vars = c("time"))
    
    # Plot number of participants at each time.
    plot_600  <- ggplot(df_plot,aes(time, value)) + 
      geom_line(aes( colour = variable))  +
      ylab("Number of participants") +
      scale_color_nejm() + 
      theme_bw()
    
    df_age1500 <- data.frame( time = z$time, ob = ob_data_t[,1501], sim = z$`1500`) 
    df_plot <- reshape2::melt(df_age1500, id.vars = c("time"))
    
    # Plot number of participants at each time.
    plot_1500  <- ggplot(df_plot,aes(time, value)) + 
      geom_line(aes( colour = variable))  +
      ylab("Number of participants") +
      scale_color_nejm() + 
      theme_bw()
    
    df_tab <- data.frame(variable = c("LL", "gamma1","gamma2","gamma3"),
                         value = mat_sort[it,])
    table.p <- ggtexttable(df_tab)
    
    plot <- ggarrange(table.p ,plot_1,plot_31,
                      plot_300,plot_600, plot_1500,
               common.legend = TRUE, legend="bottom")
    path <- paste0("~/MAD_MODEL/SUR_MODEL/Code/PARAM_ESTIMATION/MH/plot_",it, ".png")
    ggsave(path, plot= last_plot(), device = "png")
    
    print(paste("it:",it))
    plot
}

