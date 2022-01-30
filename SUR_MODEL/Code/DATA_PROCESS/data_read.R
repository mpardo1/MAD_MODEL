rm(list = ls())
library(easypackages)
libraries("ggplot2","tidyverse","ggstatsplot","deSolve")

# Mosquito reports:
Path = "~/MAD_MODEL/SUR_MODEL/data/a000_mosquito_alert_spatio_temporal_data_D_mod_df.Rds"
df <- read_rds(Path) %>% filter(presence=TRUE)

# Death rate:
Path = "~/MAD_MODEL/SUR_MODEL/data/user_participation_mortality_table.Rds"
df_mort <- readRDS(Path)

# Number of P_i:
Path_age = "~/MAD_MODEL/SUR_MODEL/data/ages_days.csv"
df_age <- read_csv(Path_age)

# Propensity probability:
Path_prop = "~/MAD_MODEL/SUR_MODEL/data/propensity_predictions.csv"
df_prop <- read_csv(Path_prop)

# Number of P_i in BCN:
Path_age = "~/MAD_MODEL/SUR_MODEL/data/ages_days_bcn.csv"
df_age_bcn <- read_csv(Path_age)


Path_life_t <- "~/MAD_MODEL/SUR_MODEL/data/participation_life_table.csv"
df_lt <- read_csv(Path_life_t)


