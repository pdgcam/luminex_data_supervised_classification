
library(ggplot2)
library(cowplot)
library(RColorBrewer)
library(matrixStats)
library(stringr)
library(data.table)
library(dplyr)
library(scales)
library(purrr)
library(tidyr)
library(stringr)
library(sf)
library(rnaturalearth)
library(rnaturalearthdata)
library(dplyr)
library(here)
library(terra)
library(exactextractr)
library(raster)
library(readxl)
library(knitr)


#--- source functions
source(here('/Users/ap2488/Documents/GitHub/updates_uminex_data_supervised_classification/Functions.R'))


# --- Import data (clean ie removed subclincal infections)
clean_HI_ratio_df <- read.csv("/Users/ap2488/Desktop/supervised_learning_flavi/FinalLuminexClassification/clean_HI_ratio_df.csv")
colnames(clean_HI_ratio_df)

# get ratios of post / pre infection 
# --- Post / pre ratio + Cross-sectional data 
all_preprocessed_dfs <- process_luminex_data(clean_HI_ratio_df, patient_pcr_mapping, pre_threshold = -1)

# ---- save pre-processed dfs ----
saveRDS(all_processed_dfs_with_chik, here("luminex_processed_data_with_chik.rds"))



# Ratio plot 