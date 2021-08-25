ex_frame <- data.frame(iterations = seq(1,10,1), 
                      ex_time = c(2.211089,2.334053, 2.75596,
                                  3.115905, 3.488181,3.647778,
                                  4.165789, 4.400713,4.926383,
                                  5.020474))
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

predict(fit, data.frame(iterations = 100000))/86400


