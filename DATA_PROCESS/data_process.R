rm(list = ls())
library(easypackages)
libraries("ggplot2","tidyverse","ggstatsplot","deSolve")

# Mosquito reports:
Path = "/Users/celsaaraujobarja/MAD_MODEL/DATA_PROCESS/Data/user_participation_mortality_table.Rds"
df <- read_rds("/Users/celsaaraujobarja/Documents/PHD/2021/JohnData/a000_mosquito_alert_spatio_temporal_data_D_mod_df.Rds") %>% filter(presence=TRUE)

# Death rate:
Path = "/Users/celsaaraujobarja/MAD_MODEL/DATA_PROCESS/Data/user_participation_mortality_table.Rds"
df_mort <- readRDS(Path)

# Number of P_i:
Path_age = "/Users/celsaaraujobarja/MAD_MODEL/DATA_PROCESS/Data/ages_days.csv"
df_age <- read_csv(Path_age)

# Propensity probability:
Path_prop = "/Users/celsaaraujobarja/MAD_MODEL/DATA_PROCESS/Data/propensity_predictions.csv"
df_prop <- read_csv(Path_prop)

# Number of P_i in BCN:
Path_age = "/Users/celsaaraujobarja/MAD_MODEL/DATA_PROCESS/Data/ages_days_bcn.csv"
df_age_bcn <- read_csv(Path_age)


Path_life_t <- "/Users/celsaaraujobarja/MAD_MODEL/DATA_PROCESS/Data/participation_life_table.csv"
df_lt <- read_csv(Path_life_t)
