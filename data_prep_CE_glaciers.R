############################
# Data preparation for analysis of the data from a DCE on Chinese public's preferences for glacier functions/ecosystem services
# n=3000, conducted online in July 2024
# study authors: Can Zhang, Bartosz Bartkowski
# code author: Bartosz Bartkowski
# contact: bartosz.bartkowski@ufz.de
# This script can be found at https://github.com/BartoszBartk/glacier-ce
############################

# packages
require(here)
require(tidyverse)

# load design dataset
design.raw <- read.csv(here("Design/design_for_analysis.csv"), sep = ";", header = T)

# load responses dataset
ce.data <- read.csv(here("Data/Glacier_CE_raw_data_20240923.csv"), sep = ";", header = T)

# recode attribute levels to categorical
# climate
design.raw$alt1 <- ifelse(design.raw$attribute == "climate" & design.raw$alt1 == 0, "1.8_degrees",
                          ifelse(design.raw$attribute == "climate" & design.raw$alt1 == 1, "1.5_degrees",
                                 ifelse(design.raw$attribute == "climate" & design.raw$alt1 == 2, "1.2_degrees",
                                        design.raw$alt1)))
design.raw$alt2 <- ifelse(design.raw$attribute == "climate" & design.raw$alt2 == 0, "1.8_degrees",
                          ifelse(design.raw$attribute == "climate" & design.raw$alt2 == 1, "1.5_degrees",
                                 ifelse(design.raw$attribute == "climate" & design.raw$alt2 == 2, "1.2_degrees",
                                        design.raw$alt2)))
# flood
design.raw$alt1 <- ifelse(design.raw$attribute == "flood" & design.raw$alt1 == 0, "frequent",
                          ifelse(design.raw$attribute == "flood" & design.raw$alt1 == 1, "10_%",
                                 ifelse(design.raw$attribute == "flood" & design.raw$alt1 == 2, "20_%",
                                        design.raw$alt1)))
design.raw$alt2 <- ifelse(design.raw$attribute == "flood" & design.raw$alt2 == 0, "frequent",
                          ifelse(design.raw$attribute == "flood" & design.raw$alt2 == 1, "10_%",
                                 ifelse(design.raw$attribute == "flood" & design.raw$alt2 == 2, "20_%",
                                        design.raw$alt2)))
# drought
design.raw$alt1 <- ifelse(design.raw$attribute == "drought" & design.raw$alt1 == 0, "frequent",
                          ifelse(design.raw$attribute == "drought" & design.raw$alt1 == 1, "10_%",
                                 ifelse(design.raw$attribute == "drought" & design.raw$alt1 == 2, "20_%",
                                        design.raw$alt1)))
design.raw$alt2 <- ifelse(design.raw$attribute == "drought" & design.raw$alt2 == 0, "frequent",
                          ifelse(design.raw$attribute == "drought" & design.raw$alt2 == 1, "10_%",
                                 ifelse(design.raw$attribute == "drought" & design.raw$alt2 == 2, "20_%",
                                        design.raw$alt2)))
# water
design.raw$alt1 <- ifelse(design.raw$attribute == "freshwater" & design.raw$alt1 == 0, "30_mil",
                          ifelse(design.raw$attribute == "freshwater" & design.raw$alt1 == 1, "35_mil",
                                 ifelse(design.raw$attribute == "freshwater" & design.raw$alt1 == 2, "40_mil",
                                        design.raw$alt1)))
design.raw$alt2 <- ifelse(design.raw$attribute == "freshwater" & design.raw$alt2 == 0, "30_mil",
                          ifelse(design.raw$attribute == "freshwater" & design.raw$alt2 == 1, "35_mil",
                                 ifelse(design.raw$attribute == "freshwater" & design.raw$alt2 == 2, "40_mil",
                                        design.raw$alt2)))
# habitat
design.raw$alt1 <- ifelse(design.raw$attribute == "habitat" & design.raw$alt1 == 0, "high_risk",
                          ifelse(design.raw$attribute == "habitat" & design.raw$alt1 == 1, "half_risk",
                                 design.raw$alt1))
design.raw$alt2 <- ifelse(design.raw$attribute == "habitat" & design.raw$alt2 == 0, "high_risk",
                          ifelse(design.raw$attribute == "habitat" & design.raw$alt2 == 1, "half_risk",
                                 design.raw$alt2))

# same for sq
design.raw$sq <- ifelse(design.raw$attribute == "climate", "1.8_degrees",
                        ifelse(design.raw$attribute == "flood", "frequent",
                               ifelse(design.raw$attribute == "drought", "frequent",
                                      ifelse(design.raw$attribute == "freshwater", "30_mil",
                                             ifelse(design.raw$attribute == "habitat", "high_risk",
                                                    design.raw$sq)))))

# reshape design dataset to have one row per choice set
design <- reshape(design.raw, 
                  direction = "wide",
                  idvar = "choice.set",
                  timevar = "attribute",
                  v.names = c("alt1", "alt2", "sq"))

# reshape responses dataset to have one row per choice set
ce_data <- reshape(ce.data,
                   direction = "long",
                   varying = c("q1", "q2", "q3", "q4", "q5", "q6", "q7", "q8", "q9", "q10"),
                   v.names = "RES",
                   timevar = "choice.set",
                   idvar = "ID")

# merge both datasets
ce_data_full <- merge(ce_data, design, by = "choice.set")

# reorder by ID
ce_data_full <- ce_data_full %>% arrange(ID)

# save full dataset
write.csv(ce_data_full, here("Data/ce_main_dataset.csv"))
