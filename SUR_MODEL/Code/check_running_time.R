rm(list = ls())
library("parallel")
library("tidyverse")
library("ggplot2")

ex_frame <- data.frame(iterations = c(seq(1,10,1),100,10000), 
                      ex_time = c(2.211089,2.334053, 2.75596,
                                  3.115905, 3.488181,3.647778,
                                  4.165789, 4.400713,4.926383,
                                  5.020474, 49.82251,8667.972))

# Two chains:
ex_frame <- data.frame(iterations = c(1,10,50,100,500), 
                       ex_time = c(2.211089,8.057199, 34.4963,
                                   67.76064,349.3564))

ex_frame <- data.frame(iterations = c(1,10,100,200), 
                       ex_time = c(0.228083833333, 0.79517883333, 6.350967, 12.52312))

ex_frame <- data.frame(iterations = c(1,10,100,150), 
                       ex_time = c(0,45304466667, 1.607306, 13.35163, 12.52312))

ggplot(ex_frame) +
 geom_line(aes(iterations, ex_time))

fit <- lm( ex_time ~ iterations, data = ex_frame)
ggplot(data=ex_frame, aes(fit$residuals)) +
  geom_histogram(binwidth = 1, color = "black", fill = "purple4") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(),
        axis.line.y=element_line()) +
  ggtitle("Histogram for Model Residuals")

ggplot(data = ex_frame, aes(x = iterations, y = ex_time)) +
  geom_point() +
  stat_smooth(method = "lm", col = "dodgerblue3") +
  theme(panel.background = element_rect(fill = "white"),
        axis.line.x=element_line(),
        axis.line.y=element_line()) +
  ggtitle("Linear Model Fitted to Data")

predict(fit, data.frame(iterations = 10000))/3600


