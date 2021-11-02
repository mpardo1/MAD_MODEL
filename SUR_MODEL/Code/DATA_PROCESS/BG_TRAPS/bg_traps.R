rm(list = ls())
library(gdata) 
library(ggplot2)
library(numbers)
library(tidyverse)
library(data.table)
library(multiplex)
library(tidyverse)
library("readxl")
library(reshape)
library(viridis)
library(stats)

# Upload data from bgtraps.
Path = '/home/marta/Documentos/PHD/2021/SUR_Model/Code/data/bcn_bgcount_time_profile.csv'
bg_traps= read.csv(Path)
ggplot(bg_traps, aes(x = date, y = value)) +
  geom_point()
density = 1664182
bg_traps$prop <- bg_traps$value/density
ggplot(bg_traps, aes(x = date, y = prop)) +
  geom_line()
