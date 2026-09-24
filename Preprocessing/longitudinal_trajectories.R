
library(here)

#--- source functions
source(here('/Users/ap2488/Documents/GitHub/luminex_data_supervised_classification/Functions.R'))

# --- Read data 
# validation and random datasets (cebu) - validation subset = PCR confirmed cases
validation_subset <- read.csv(here("Data/MIA_DataBaseOut_ValidationSet.csv"))
random_subset <- read.csv(here("Data/MIA_DataBaseOut_RandomSubset.csv"))
cebu_mutiple_antigens <- read_excel(here("Data/db_philippines_IgG_IgA_IgM_avidity.xlsx"))
cpc_gps <- read.csv(here("Data/CPC_GPS.csv"))

head(cebu_mutiple_antigens)

# align patient IDs / PCR cols across datasets
cebu_mutiple_antigens$id_patient <- gsub("_", "-", cebu_mutiple_antigens$id_patient)
length(intersect(validation_subset$ids, cebu_mutiple_antigens$id_patient)) #39 samples intersect


# pivot to get RAU as main data
cebu_pivot <- cebu_mutiple_antigens %>%
  pivot_wider(
    names_from = antigen,
    values_from = RAU
  )

# --- Add col: days_since_infection 
cebu_pivot_days_since_inf <- cebu_pivot %>%
  group_by(id_patient) %>%
  arrange(date_sample) %>%
  mutate(
    n_samples = n(),
    first_positive_date = first(date_sample[PCR %in% c("DENV1", "DENV2", "DENV3", "DENV4", "ZIKV", "CHIKV")],
                                 default = as.Date(NA)),
    first_negative_date = first(date_sample[PCR %in% "negative"],
                                 default = as.Date(NA)),
    infection_date = if_else(!is.na(first_positive_date), first_positive_date, first_negative_date),
    days_since_infection = as.numeric(difftime(date_sample, infection_date, units = "days"))
  ) %>%
  ungroup() %>%
  dplyr::select(-first_positive_date, -first_negative_date)


colnames(cebu_pivot_days_since_inf)
head(cebu_pivot_days_since_inf$DENV3_NS1)

# look at NS1 trajectory for DENV1, DENV2, DENV3, DENV4, ZIKV, YF, WNV, JEV (IgG)
ns1_cols <- c("DENV1_NS1","DENV2_NS1","DENV3_NS1","DENV4_NS1",
              "JEV_NS1","WNV_NS1","YFV_NS1","ZIKV_NS1","ZIKVSU_NS1")
ns1_traj_df <- cebu_pivot_days_since_inf %>%
  filter(isotype == "IgG") %>%
  dplyr::select(id_patient, PCR, days_since_infection, n_samples, all_of(ns1_cols)) %>%
  pivot_longer(all_of(ns1_cols), names_to = "antigen", values_to = "RAU") %>%
  filter(!is.na(RAU), !is.na(days_since_infection)) %>%
  mutate(
    virus  = sub("_NS1$", "", antigen),
    family = sub("^ZIKV.*$", "ZIKV", virus),
    family = factor(family,
      levels = c("DENV1","DENV2","DENV3","DENV4","JEV","WNV","YFV","ZIKV"))
  )
n_ids <- n_distinct(ns1_traj_df$id_patient)

fam_cols <- c(
  DENV1 = "#BDD7E7",   # light blue
  DENV2 = "#6BAED6",
  DENV3 = "#3182BD",
  DENV4 = "#08519C",   # dark blue
  JEV   = "#E6550D",   # orange
  WNV   = "#31A354",   # green
  YFV   = "#DAA520",   # goldenrod
  ZIKV  = "#B0179C"    # magenta
)


ns1_traj_df$RAU <- log10(ns1_traj_df$RAU)
quartz()
ns1_traj <- ggplot(ns1_traj_df,
                   aes(days_since_infection, RAU,
                       colour = family, group = antigen)) +
  geom_line(linewidth = 0.4) +
  geom_point(size = 0.7) +
  facet_wrap(~ id_patient, ncol = 8) +
  scale_y_log10() +
  scale_colour_manual(values = fam_cols, name = "Virus") +
  labs(x = "Days since infection", y = "NS1 IgG (RAU)") +
  theme_minimal(base_size = 8) +
  theme(panel.grid.minor = element_blank(),
        strip.text = element_text(face = "bold", size = 7),
        legend.position = "bottom")
print(ns1_traj)
