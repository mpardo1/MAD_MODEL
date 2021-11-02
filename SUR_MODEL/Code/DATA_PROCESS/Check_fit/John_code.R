rm(list = ls())
library(gdata) 
library(ggplot2)
library(numbers)
library(tidyverse)
library(data.table)


agepyr = D %>% filter(!is.na(age), !is.na(gender)) %>% 
  group_by(age, gender, wave) %>% summarise(n = n()) %>% ungroup() 
  %>% mutate(wave = paste0("Wave ", wave), n = case_when(gender=="female"~n, gender=="male"~-n))

p = ggplot(agepyr, aes(x = age, y = n, fill = gender, text=paste('Age: ', age, '<br>Gender: ', gender, '<br>Number: ', abs(n)))) 
+ geom_bar(data = subset(agepyr, gender == "female"), stat = "identity") 
+ geom_bar(data = subset(agepyr, gender == "male"), stat = "identity") 
+  coord_flip() 
+ scale_y_continuous(breaks = seq(-80, 80, 20), labels = as.character(c(seq(80,0, -20), seq(20, 80, 20))), limits = c(-80, 80)) 
+ ylab("number of respondents") + facet_grid(~wave) 
+ theme(legend.title = element_blank())

gp = ggplotly(p, tooltip = "text") %>% layout(legend = list(title = list(text="gender")
    , orientation = "v", yanchor = "center", y = 0.5)) %>% config(displaylogo = FALSE)

frameWidget(gp)

library(lubridate)

x <- c("02/01/2000", "20/02/2000", "12/12/2000", "13/01/2001")
date <- dmy(x)

days <- yday(date) - 1 # so Jan 1 = day 0 
total_days <- cumsum(days)