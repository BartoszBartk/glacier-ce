############################
# Analysis of the data from a DCE on Chinese public's preferences for glacier functions/ecosystem services
# n=3000, conducted online in July 2024
# study authors: Can Zhang, Bartosz Bartkowski
# code author: Bartosz Bartkowski
# contact: bartosz.bartkowski@ufz.de
# This script can be found at https://github.com/BartoszBartk/glacier-ce
############################

############## Outline ###############
# 1. Packages and data
# 2. Multinomial logit
# 3. Simple mixed logit
# 4. Mixed logit with regional interactions
# 5. Simple mixed logit in WTP space
# 6. Mixed logit with all attribute interactions in WTP space
# 7. Mixed logit with SQ interactions in WTP space
# 8. Mixed logit with all relevant interactions in preference space
# 9. Multinomial logit with all relevant interactions in preference space
# 10. Mixed logit with regional interactions in WTP space
# 11. [test] Multinomial logit with continuous variables
# 12. Save outputs as DOCX table(s)
######################################

######################
# 1. Packages and data
######################

# optional: clear environment
rm(list=ls())

# load packages
require(apollo)
require(here)
require(choiceTools)
require(texreg)

# read in dataset already formatted for analysis
# the dataset used here has already been checked for protest votes (none found)
database <- read.csv(here("Data/ce_main_dataset.csv"), header = T)

# add some modifications to database to simplify further analysis
# make experience variables binary
database$exp_drought <- ifelse(database$exp_drought > 0, 1, 0)
database$exp_flood <- ifelse(database$exp_flood > 0, 1, 0)
# only use highest ranked means of protecting glaciers
database$way_climate.mitigation <- ifelse(database$way_climate.mitigation != 1 | is.na(database$way_climate.mitigation), 0, 1)
database$way_national.park <- ifelse(database$way_national.park != 1 | is.na(database$way_national.park), 0, 1)
database$way_pollution.control <- ifelse(database$way_pollution.control != 1 | is.na(database$way_pollution.control), 0, 1)
# make income variable pseudo-continuous
database$income <- ifelse(database$income == "10000_below", 8000,
                          ifelse(database$income == "10000_20000", 15000,
                                 ifelse(database$income == "20000_40000", 30000,
                                        ifelse(database$income == "40000_60000", 50000,
                                               ifelse(database$income == "60000_80000", 90000,
                                                      ifelse(database$income == "80000_100000", 90000,
                                                             ifelse(database$income == "100000_120000", 110000,
                                                                    ifelse(database$income == "120000_140000", 130000,
                                                                           ifelse(database$income == "140000_160000", 150000,
                                                                                  ifelse(database$income == "160000_200000", 180000,
                                                                                         ifelse(database$income == "200000_above", 250000, database$income)))))))))))
database$income <- as.numeric(database$income)
# make education variable binary (higher education vs. below)
database$edu <- ifelse(database$edu == "higher_education", 1, 0)
database$edu <- as.numeric(database$edu)

# initialize apollo package
apollo_initialise()

######################
# 2. Multinomial logit
######################

## implement a simple MNL
# set required parameters of project
apollo_control <- list(
  modelName="glaciers_CE_mnl",
  modelDescr="Discrete choice experiment into preferences for glacier ES in China",
  indivID="ID"
)

# define coefficients to be estimated and starting values for each
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = 0,
                 b_climate_1.8 = 0,
                 b_climate_1.5 = 0,
                 b_climate_1.2 = 0,
                 b_flood_freq = 0,
                 b_flood_10 = 0,
                 b_flood_20 = 0,
                 b_drought_freq = 0,
                 b_drought_10 = 0,
                 b_drought_20 = 0,
                 b_water_30 = 0,
                 b_water_35 = 0,
                 b_water_40 = 0,
                 b_habitat_high = 0,
                 b_habitat_half = 0,
                 b_cost = 0)

#set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high")

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  V = list()
  # define utilities for each alternative (brackets required because of the breaks)
  V[['alt1']] = (asc_1 + 
    b_climate_1.8 * (alt1.climate == "1.8_degrees") + b_climate_1.5 * (alt1.climate == "1.5_degrees") + b_climate_1.2 * (alt1.climate == "1.2_degrees") +
    b_flood_freq * (alt1.flood == "frequent") + b_flood_10 * (alt1.flood == "10_%") + b_flood_20 * (alt1.flood == "20_%") +
    b_drought_freq * (alt1.drought == "frequent") + b_drought_10 * (alt1.drought == "10_%") + b_drought_20 * (alt1.drought == "20_%") + 
    b_water_30 * (alt1.freshwater == "30_mil") + b_water_35 * (alt1.freshwater == "35_mil") + b_water_40 * (alt1.freshwater == "40_mil") +
    b_habitat_high * (alt1.habitat == "high_risk") + b_habitat_half * (alt1.habitat == "half_risk") +
    b_cost * alt1.cost)
  V[['alt2']] = (asc_2 + 
    b_climate_1.8 * (alt2.climate == "1.8_degrees") + b_climate_1.5 * (alt2.climate == "1.5_degrees") + b_climate_1.2 * (alt2.climate == "1.2_degrees") +
    b_flood_freq * (alt2.flood == "frequent") + b_flood_10 * (alt2.flood == "10_%") + b_flood_20 * (alt2.flood == "20_%") +
    b_drought_freq * (alt2.drought == "frequent") + b_drought_10 * (alt2.drought == "10_%") + b_drought_20 * (alt2.drought == "20_%") + 
    b_water_30 * (alt2.freshwater == "30_mil") + b_water_35 * (alt2.freshwater == "35_mil") + b_water_40 * (alt2.freshwater == "40_mil") +
    b_habitat_high * (alt2.habitat == "high_risk") + b_habitat_half * (alt2.habitat == "half_risk") +
    b_cost * alt2.cost)
  V[['alt3']] = (asc_3 + 
    b_climate_1.8 * (sq.climate == "1.8_degrees") +# b_climate_1.5 * (sq.climate == "1.5_degrees") + b_climate_1.2 * (sq.climate == "1.2_degrees") +
    b_flood_freq * (sq.flood == "frequent") +# b_flood_10 * (sq.flood == "10_%") + b_flood_20 * (sq.flood == "20_%") +
    b_drought_freq * (sq.drought == "frequent") +# b_drought_10 * (sq.drought == "10_%") + b_drought_20 * (sq.drought == "20_%") + 
    b_water_30 * (sq.freshwater == "30_mil") +# b_water_35 * (sq.freshwater == "35_mil") + b_water_40 * (sq.freshwater == "40_mil") +
    b_habitat_high * (sq.habitat == "high_risk") +# b_habitat_half * (sq.habitat == "half_risk") +
    b_cost * sq.cost)
    
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

#estimate model
mnl <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show outputs
apollo_modelOutput(mnl,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mnl,
                  saveOutput_settings=list(printPVal=1))

#######################
# 3. Simple mixed logit
#######################

## implement a simple mixed logit model to allow for inter-individual preference heterogeneity
apollo_control = list(
  modelName = "glaciers_CE_mxl",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  mixing = TRUE,
  nCores = 7,
  seed = 2104
)

# define coefficients to be estimated and starting values for each
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = 0,
                 b_climate_1.8 = 0,
                 m_b_climate_1.5 = 0.01,
                 m_b_climate_1.2 = 0.01,
                 b_flood_freq = 0,
                 m_b_flood_10 = 0.01,
                 m_b_flood_20 = 0.01,
                 b_drought_freq = 0,
                 m_b_drought_10 = 0.01,
                 m_b_drought_20 = 0.01,
                 b_water_30 = 0,
                 m_b_water_35 = 0.01,
                 m_b_water_40 = 0.01,
                 b_habitat_high = 0,
                 m_b_habitat_half = 0.01,
                 m_log_b_cost = -1,
                 #sigma_b_climate_1.8 = 0,
                 sigma_b_climate_1.5 = 0.01,
                 sigma_b_climate_1.2 = 0.01,
                 #sigma_b_flood_freq = 0,
                 sigma_b_flood_10 = 0.01,
                 sigma_b_flood_20 = 0.01,
                 #sigma_b_drought_freq = 0,
                 sigma_b_drought_10 = 0.01,
                 sigma_b_drought_20 = 0.01,
                 #sigma_b_water_30 = 0,
                 sigma_b_water_35 = 0.01,
                 sigma_b_water_40 = 0.01,
                 #sigma_b_habitat_high = 0,
                 sigma_b_habitat_half = 0.01,
                 sigma_log_b_cost = 1)

#set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high"
)

#define parameters for the simulation
apollo_draws = list(
  interDrawsType = "sobol",
  interNDraws = 1000,
  interNormDraws = c("draws_cost_inter",
                     #"draws_climate18_inter",
                     "draws_climate15_inter",
                     "draws_climate12_inter",
                     #"draws_floodfreq_inter",
                     "draws_flood10_inter",
                     "draws_flood20_inter",
                     #"draws_droughtfreq_inter",
                     "draws_drought10_inter",
                     "draws_drought20_inter",
                     #"draws_water30_inter",
                     "draws_water35_inter",
                     "draws_water40_inter",
                     #"draws_habitathi_inter",
                     "draws_habitatha_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_climate_1.5"]] = m_b_climate_1.5 + sigma_b_climate_1.5 * draws_climate15_inter
  randcoeff[["b_climate_1.2"]] = m_b_climate_1.2 + sigma_b_climate_1.2 * draws_climate12_inter
  randcoeff[["b_flood_10"]] = m_b_flood_10 + sigma_b_flood_10 * draws_flood10_inter
  randcoeff[["b_flood_20"]] = m_b_flood_20 + sigma_b_flood_20 * draws_flood20_inter
  randcoeff[["b_drought_10"]] = m_b_drought_10 + sigma_b_drought_10 * draws_drought10_inter
  randcoeff[["b_drought_20"]] = m_b_drought_20 + sigma_b_drought_20 * draws_drought20_inter
  randcoeff[["b_water_35"]] = m_b_water_35 + sigma_b_water_35 * draws_water35_inter
  randcoeff[["b_water_40"]] = m_b_water_40 + sigma_b_water_40 * draws_water40_inter
  randcoeff[["b_habitat_half"]] = m_b_habitat_half + sigma_b_habitat_half * draws_habitatha_inter
  randcoeff[["b_cost"]] = -exp(m_log_b_cost + sigma_log_b_cost * draws_cost_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details, see apollo documentation)
apollo_probabilities = function(apollo_beta,apollo_inputs,functionality="estimate"){
  apollo_attach(apollo_beta,apollo_inputs)
  on.exit(apollo_detach(apollo_beta,apollo_inputs))
  P = list()
  V = list()
  V[['alt1']] = (asc_1 + 
                   b_climate_1.8 * (alt1.climate == "1.8_degrees") + b_climate_1.5 * (alt1.climate == "1.5_degrees") + b_climate_1.2 * (alt1.climate == "1.2_degrees") +
                   b_flood_freq * (alt1.flood == "frequent") + b_flood_10 * (alt1.flood == "10_%") + b_flood_20 * (alt1.flood == "20_%") +
                   b_drought_freq * (alt1.drought == "frequent") + b_drought_10 * (alt1.drought == "10_%") + b_drought_20 * (alt1.drought == "20_%") + 
                   b_water_30 * (alt1.freshwater == "30_mil") + b_water_35 * (alt1.freshwater == "35_mil") + b_water_40 * (alt1.freshwater == "40_mil") +
                   b_habitat_high * (alt1.habitat == "high_risk") + b_habitat_half * (alt1.habitat == "half_risk") +
                   b_cost * alt1.cost)
  V[['alt2']] = (asc_2 + 
                   b_climate_1.8 * (alt2.climate == "1.8_degrees") + b_climate_1.5 * (alt2.climate == "1.5_degrees") + b_climate_1.2 * (alt2.climate == "1.2_degrees") +
                   b_flood_freq * (alt2.flood == "frequent") + b_flood_10 * (alt2.flood == "10_%") + b_flood_20 * (alt2.flood == "20_%") +
                   b_drought_freq * (alt2.drought == "frequent") + b_drought_10 * (alt2.drought == "10_%") + b_drought_20 * (alt2.drought == "20_%") + 
                   b_water_30 * (alt2.freshwater == "30_mil") + b_water_35 * (alt2.freshwater == "35_mil") + b_water_40 * (alt2.freshwater == "40_mil") +
                   b_habitat_high * (alt2.habitat == "high_risk") + b_habitat_half * (alt2.habitat == "half_risk") +
                   b_cost * alt2.cost)
  V[['alt3']] = (asc_3 + 
                   b_climate_1.8 * (sq.climate == "1.8_degrees") +# b_climate_1.5 * (sq.climate == "1.5_degrees") + b_climate_1.2 * (sq.climate == "1.2_degrees") +
                   b_flood_freq * (sq.flood == "frequent") +# b_flood_10 * (sq.flood == "10_%") + b_flood_20 * (sq.flood == "20_%") +
                   b_drought_freq * (sq.drought == "frequent") +# b_drought_10 * (sq.drought == "10_%") + b_drought_20 * (sq.drought == "20_%") + 
                   b_water_30 * (sq.freshwater == "30_mil") +# b_water_35 * (sq.freshwater == "35_mil") + b_water_40 * (sq.freshwater == "40_mil") +
                   b_habitat_high * (sq.habitat == "high_risk") +# b_habitat_half * (sq.habitat == "half_risk") +
                   b_cost * sq.cost)
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings,functionality)
  P = apollo_panelProd(P,apollo_inputs,functionality)
  P = apollo_avgInterDraws(P,apollo_inputs,functionality)
  P = apollo_prepareProb(P,apollo_inputs,functionality)
  return(P)
}

#estimate model
mxl <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show output
apollo_modelOutput(mxl,
                   modelOutput_settings=list(printPVal=1))
#write output
apollo_saveOutput(mxl,
                  saveOutput_settings=list(printPVal=1))

## Calculate WTP with help of Delta method
mxl <- apollo_loadModel("glaciers_CE_mxl") # the working directory needs to be the one where the model is stored for this to work

# optionally: calculate mean and SD of the lognormal distribution of the price parameter
deltaMethod_settings = list(operation = "lognormal",
                            parName1 = "m_log_b_cost",
                            parName2 = "sigma_log_b_cost")
apollo_deltaMethod(mxl, deltaMethod_settings)

# calculate WTP mean (with formula for moments of the price parameter distribution written out)
deltaMethod_settings = list(expression = c(wtp_climate15 = "m_b_climate_1.5/exp(m_log_b_cost + (sigma_log_b_cost^2)/2)",
                                           wtp_climate12 = "m_b_climate_1.2/exp(m_log_b_cost + (sigma_log_b_cost^2)/2)",
                                           wtp_flood10 = "m_b_flood_10/exp(m_log_b_cost + (sigma_log_b_cost^2)/2)",
                                           wtp_flood20 = "m_b_flood_20/exp(m_log_b_cost + (sigma_log_b_cost^2)/2)",
                                           wtp_drought10 = "m_b_drought_10/exp(m_log_b_cost + (sigma_log_b_cost^2)/2)",
                                           wtp_drought20 = "m_b_drought_20/exp(m_log_b_cost + (sigma_log_b_cost^2)/2)",
                                           wtp_water35 = "m_b_water_35/exp(m_log_b_cost + (sigma_log_b_cost^2)/2)",
                                           wtp_water40 = "m_b_water_40/exp(m_log_b_cost + (sigma_log_b_cost^2)/2)",
                                           wtp_habhalf = "m_b_habitat_half/exp(m_log_b_cost + (sigma_log_b_cost^2)/2)"
                                           ))
apollo_deltaMethod(mxl, deltaMethod_settings)

######################
# 4. Mixed logit with regional interactions
######################

apollo_control <- list(
  modelName="glaciers_CE_mxl_reg",
  modelDescr="Discrete choice experiment into preferences for glacier ES in China",
  indivID="ID",
  mixing = TRUE,
  nCores = 7,
  seed = 2104
)

# define coefficients to be estimated and starting values for each
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = 0,
                 b_climate_1.8 = 0,
                 m_b_climate_1.5 = 0.01,
                 m_b_climate_1.2 = 0.01,
                 b_flood_freq = 0,
                 m_b_flood_10 = 0.01,
                 m_b_flood_20 = 0.01,
                 b_drought_freq = 0,
                 m_b_drought_10 = 0.01,
                 m_b_drought_20 = 0.01,
                 b_water_30 = 0,
                 m_b_water_35 = 0.01,
                 m_b_water_40 = 0.01,
                 b_habitat_high = 0,
                 m_b_habitat_half = 0.01,
                 m_log_b_cost = -1,
                 #sigma_b_climate_1.8 = 0,
                 sigma_b_climate_1.5 = 0.01,
                 sigma_b_climate_1.2 = 0.01,
                 #sigma_b_flood_freq = 0,
                 sigma_b_flood_10 = 0.01,
                 sigma_b_flood_20 = 0.01,
                 #sigma_b_drought_freq = 0,
                 sigma_b_drought_10 = 0.01,
                 sigma_b_drought_20 = 0.01,
                 #sigma_b_water_30 = 0,
                 sigma_b_water_35 = 0.01,
                 sigma_b_water_40 = 0.01,
                 #sigma_b_habitat_high = 0,
                 sigma_b_habitat_half = 0.01,
                 sigma_log_b_cost = 1,
                 b_basinds_climate15 = 0,
                 b_basinds_climate12 = 0,
                 b_basinds_flood10 = 0,
                 b_basinds_flood20 = 0,
                 b_basinds_drought10 = 0,
                 b_basinds_drought20 = 0,
                 b_basinds_water35 = 0,
                 b_basinds_water40 = 0,
                 b_basinds_habitathalf = 0,
                 b_basinwt_climate15 = 0,
                 b_basinwt_climate12 = 0,
                 b_basinwt_flood10 = 0,
                 b_basinwt_flood20 = 0,
                 b_basinwt_drought10 = 0,
                 b_basinwt_drought20 = 0,
                 b_basinwt_water35 = 0,
                 b_basinwt_water40 = 0,
                 b_basinwt_habitathalf = 0)

#set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high"
                  )

#define parameters for the simulation
apollo_draws = list(interDrawsType = "sobol",
                    interNDraws = 1000,
                    interNormDraws = c("draws_cost_inter",
                                      #"draws_climate18_inter",
                                      "draws_climate15_inter",
                                      "draws_climate12_inter",
                                      #"draws_floodfreq_inter",
                                      "draws_flood10_inter",
                                      "draws_flood20_inter",
                                      #"draws_droughtfreq_inter",
                                      "draws_drought10_inter",
                                      "draws_drought20_inter",
                                      #"draws_water30_inter",
                                      "draws_water35_inter",
                                      "draws_water40_inter",
                                      #"draws_habitathi_inter",
                                      "draws_habitatha_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_climate_1.5"]] = m_b_climate_1.5 + sigma_b_climate_1.5 * draws_climate15_inter
  randcoeff[["b_climate_1.2"]] = m_b_climate_1.2 + sigma_b_climate_1.2 * draws_climate12_inter
  randcoeff[["b_flood_10"]] = m_b_flood_10 + sigma_b_flood_10 * draws_flood10_inter
  randcoeff[["b_flood_20"]] = m_b_flood_20 + sigma_b_flood_20 * draws_flood20_inter
  randcoeff[["b_drought_10"]] = m_b_drought_10 + sigma_b_drought_10 * draws_drought10_inter
  randcoeff[["b_drought_20"]] = m_b_drought_20 + sigma_b_drought_20 * draws_drought20_inter
  randcoeff[["b_water_35"]] = m_b_water_35 + sigma_b_water_35 * draws_water35_inter
  randcoeff[["b_water_40"]] = m_b_water_40 + sigma_b_water_40 * draws_water40_inter
  randcoeff[["b_habitat_half"]] = m_b_habitat_half + sigma_b_habitat_half * draws_habitatha_inter
  randcoeff[["b_cost"]] = -exp(m_log_b_cost + sigma_log_b_cost * draws_cost_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  
  # define interaction terms
  #b_climate18 = b_climate_1.8 + b_basinoc_climate18 * (basin == "oc") + b_basinds_climate18 * (basin == "ds") + b_basinwt_climate18 * (basin == "wt")
  b_climate15 = b_climate_1.5 + b_basinds_climate15 * (basin == "ds") + b_basinwt_climate15 * (basin == "wt")
  b_climate12 = b_climate_1.2 + b_basinds_climate12 * (basin == "ds") + b_basinwt_climate12 * (basin == "wt")
  #b_floodfreq = b_flood_freq + b_basinoc_floodfreq * (basin == "oc") + b_basinds_floodfreq * (basin == "ds") + b_basinwt_floodfreq * (basin == "wt")
  b_flood10 = b_flood_10 + b_basinds_flood10 * (basin == "ds") + b_basinwt_flood10 * (basin == "wt")
  b_flood20 = b_flood_20 + b_basinds_flood20 * (basin == "ds") + b_basinwt_flood20 * (basin == "wt")
  #b_droughtfreq = b_drought_freq + b_basinoc_droughtfreq * (basin == "oc") + b_basinds_droughtfreq * (basin == "ds") + b_basinwt_droughtfreq * (basin == "wt")
  b_drought10 = b_drought_10 + b_basinds_drought10 * (basin == "ds") + b_basinwt_drought10 * (basin == "wt")
  b_drought20 = b_drought_20 + b_basinds_drought20 * (basin == "ds") + b_basinwt_drought20 * (basin == "wt")
  #b_water30 = b_water_30 + b_basinoc_water30 * (basin == "oc") + b_basinds_water30 * (basin == "ds") + b_basinwt_water30 * (basin == "wt")
  b_water35 = b_water_35 + b_basinds_water35 * (basin == "ds") + b_basinwt_water35 * (basin == "wt")
  b_water40 = b_water_40 + b_basinds_water40 * (basin == "ds") + b_basinwt_water40 * (basin == "wt")
  #b_habitathigh = b_habitat_high + b_basinoc_habitathigh * (basin == "oc") + b_basinds_habitathigh * (basin == "ds") + b_basinwt_habitathigh * (basin == "wt")
  b_habitathalf = b_habitat_half + b_basinds_habitathalf * (basin == "ds") + b_basinwt_habitathalf * (basin == "wt")
  
  V = list()
  
  V[['alt1']] = (asc_1 + 
                   b_climate_1.8 * (alt1.climate == "1.8_degrees") + b_climate15 * (alt1.climate == "1.5_degrees") + b_climate12 * (alt1.climate == "1.2_degrees") +
                   b_flood_freq * (alt1.flood == "frequent") + b_flood10 * (alt1.flood == "10_%") + b_flood20 * (alt1.flood == "20_%") +
                   b_drought_freq * (alt1.drought == "frequent") + b_drought10 * (alt1.drought == "10_%") + b_drought20 * (alt1.drought == "20_%") + 
                   b_water_30 * (alt1.freshwater == "30_mil") + b_water35 * (alt1.freshwater == "35_mil") + b_water40 * (alt1.freshwater == "40_mil") +
                   b_habitat_high * (alt1.habitat == "high_risk") + b_habitathalf * (alt1.habitat == "half_risk") +
                   b_cost * alt1.cost)
  V[['alt2']] = (asc_2 + 
                   b_climate_1.8 * (alt2.climate == "1.8_degrees") + b_climate15 * (alt2.climate == "1.5_degrees") + b_climate12 * (alt2.climate == "1.2_degrees") +
                   b_flood_freq * (alt2.flood == "frequent") + b_flood10 * (alt2.flood == "10_%") + b_flood20 * (alt2.flood == "20_%") +
                   b_drought_freq * (alt2.drought == "frequent") + b_drought10 * (alt2.drought == "10_%") + b_drought20 * (alt2.drought == "20_%") + 
                   b_water_30 * (alt2.freshwater == "30_mil") + b_water35 * (alt2.freshwater == "35_mil") + b_water40 * (alt2.freshwater == "40_mil") +
                   b_habitat_high * (alt2.habitat == "high_risk") + b_habitathalf * (alt2.habitat == "half_risk") +
                   b_cost * alt2.cost)
  V[['alt3']] = (asc_3 + 
                   b_climate_1.8 * (sq.climate == "1.8_degrees") +# b_climate_1.5 * (sq.climate == "1.5_degrees") + b_climate_1.2 * (sq.climate == "1.2_degrees") +
                   b_flood_freq * (sq.flood == "frequent") +# b_flood_10 * (sq.flood == "10_%") + b_flood_20 * (sq.flood == "20_%") +
                   b_drought_freq * (sq.drought == "frequent") +# b_drought_10 * (sq.drought == "10_%") + b_drought_20 * (sq.drought == "20_%") + 
                   b_water_30 * (sq.freshwater == "30_mil") +# b_water_35 * (sq.freshwater == "35_mil") + b_water_40 * (sq.freshwater == "40_mil") +
                   b_habitat_high * (sq.habitat == "high_risk") +# b_habitat_half * (sq.habitat == "half_risk") +
                   b_cost * sq.cost)
  
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings,functionality)
  P = apollo_panelProd(P,apollo_inputs,functionality)
  P = apollo_avgInterDraws(P,apollo_inputs,functionality)
  P = apollo_prepareProb(P,apollo_inputs,functionality)
  return(P)
}

#estimate model
mxl_reg <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show outputs
apollo_modelOutput(mxl_reg,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mxl_reg,
                  saveOutput_settings=list(printPVal=1))

#######################
# 5. Simple mixed logit in WTP space
#######################

## implement a simple mixed logit model to allow for inter-individual preference heterogeneity
apollo_control = list(
  modelName = "glaciers_CE_mxl_wtp",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  mixing = TRUE,
  nCores = 7,
  seed = 2104
)

# define coefficients to be estimated and starting values for each
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = 0,
                 b_climate_1.8 = 0,
                 m_b_climate_1.5 = 0.01,
                 m_b_climate_1.2 = 0.01,
                 b_flood_freq = 0,
                 m_b_flood_10 = 0.01,
                 m_b_flood_20 = 0.01,
                 b_drought_freq = 0,
                 m_b_drought_10 = 0.01,
                 m_b_drought_20 = 0.01,
                 b_water_30 = 0,
                 m_b_water_35 = 0.01,
                 m_b_water_40 = 0.01,
                 b_habitat_high = 0,
                 m_b_habitat_half = 0.01,
                 m_log_b_cost = -1,
                 #sigma_b_climate_1.8 = 0,
                 sigma_b_climate_1.5 = 0.01,
                 sigma_b_climate_1.2 = 0.01,
                 #sigma_b_flood_freq = 0,
                 sigma_b_flood_10 = 0.01,
                 sigma_b_flood_20 = 0.01,
                 #sigma_b_drought_freq = 0,
                 sigma_b_drought_10 = 0.01,
                 sigma_b_drought_20 = 0.01,
                 #sigma_b_water_30 = 0,
                 sigma_b_water_35 = 0.01,
                 sigma_b_water_40 = 0.01,
                 #sigma_b_habitat_high = 0,
                 sigma_b_habitat_half = 0.01,
                 sigma_log_b_cost = 1)

#set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high"
)

#define parameters for the simulation
apollo_draws = list(
  interDrawsType = "sobol",
  interNDraws = 1000,
  interNormDraws = c("draws_cost_inter",
                     #"draws_climate18_inter",
                     "draws_climate15_inter",
                     "draws_climate12_inter",
                     #"draws_floodfreq_inter",
                     "draws_flood10_inter",
                     "draws_flood20_inter",
                     #"draws_droughtfreq_inter",
                     "draws_drought10_inter",
                     "draws_drought20_inter",
                     #"draws_water30_inter",
                     "draws_water35_inter",
                     "draws_water40_inter",
                     #"draws_habitathi_inter",
                     "draws_habitatha_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_climate_1.5"]] = m_b_climate_1.5 + sigma_b_climate_1.5 * draws_climate15_inter
  randcoeff[["b_climate_1.2"]] = m_b_climate_1.2 + sigma_b_climate_1.2 * draws_climate12_inter
  randcoeff[["b_flood_10"]] = m_b_flood_10 + sigma_b_flood_10 * draws_flood10_inter
  randcoeff[["b_flood_20"]] = m_b_flood_20 + sigma_b_flood_20 * draws_flood20_inter
  randcoeff[["b_drought_10"]] = m_b_drought_10 + sigma_b_drought_10 * draws_drought10_inter
  randcoeff[["b_drought_20"]] = m_b_drought_20 + sigma_b_drought_20 * draws_drought20_inter
  randcoeff[["b_water_35"]] = m_b_water_35 + sigma_b_water_35 * draws_water35_inter
  randcoeff[["b_water_40"]] = m_b_water_40 + sigma_b_water_40 * draws_water40_inter
  randcoeff[["b_habitat_half"]] = m_b_habitat_half + sigma_b_habitat_half * draws_habitatha_inter
  randcoeff[["b_cost"]] = -exp(m_log_b_cost + sigma_log_b_cost * draws_cost_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details, see apollo documentation)
apollo_probabilities = function(apollo_beta,apollo_inputs,functionality="estimate"){
  apollo_attach(apollo_beta,apollo_inputs)
  on.exit(apollo_detach(apollo_beta,apollo_inputs))
  P = list()
  V = list()
  V[['alt1']] = (asc_1 + b_cost * (
                   b_climate_1.8 * (alt1.climate == "1.8_degrees") + b_climate_1.5 * (alt1.climate == "1.5_degrees") + b_climate_1.2 * (alt1.climate == "1.2_degrees") +
                   b_flood_freq * (alt1.flood == "frequent") + b_flood_10 * (alt1.flood == "10_%") + b_flood_20 * (alt1.flood == "20_%") +
                   b_drought_freq * (alt1.drought == "frequent") + b_drought_10 * (alt1.drought == "10_%") + b_drought_20 * (alt1.drought == "20_%") + 
                   b_water_30 * (alt1.freshwater == "30_mil") + b_water_35 * (alt1.freshwater == "35_mil") + b_water_40 * (alt1.freshwater == "40_mil") +
                   b_habitat_high * (alt1.habitat == "high_risk") + b_habitat_half * (alt1.habitat == "half_risk") +
                   alt1.cost))
  V[['alt2']] = (asc_2 + b_cost * (
                   b_climate_1.8 * (alt2.climate == "1.8_degrees") + b_climate_1.5 * (alt2.climate == "1.5_degrees") + b_climate_1.2 * (alt2.climate == "1.2_degrees") +
                   b_flood_freq * (alt2.flood == "frequent") + b_flood_10 * (alt2.flood == "10_%") + b_flood_20 * (alt2.flood == "20_%") +
                   b_drought_freq * (alt2.drought == "frequent") + b_drought_10 * (alt2.drought == "10_%") + b_drought_20 * (alt2.drought == "20_%") + 
                   b_water_30 * (alt2.freshwater == "30_mil") + b_water_35 * (alt2.freshwater == "35_mil") + b_water_40 * (alt2.freshwater == "40_mil") +
                   b_habitat_high * (alt2.habitat == "high_risk") + b_habitat_half * (alt2.habitat == "half_risk") +
                   alt2.cost))
  V[['alt3']] = (asc_3 + b_cost* (
                   b_climate_1.8 * (sq.climate == "1.8_degrees") +# b_climate_1.5 * (sq.climate == "1.5_degrees") + b_climate_1.2 * (sq.climate == "1.2_degrees") +
                   b_flood_freq * (sq.flood == "frequent") +# b_flood_10 * (sq.flood == "10_%") + b_flood_20 * (sq.flood == "20_%") +
                   b_drought_freq * (sq.drought == "frequent") +# b_drought_10 * (sq.drought == "10_%") + b_drought_20 * (sq.drought == "20_%") + 
                   b_water_30 * (sq.freshwater == "30_mil") +# b_water_35 * (sq.freshwater == "35_mil") + b_water_40 * (sq.freshwater == "40_mil") +
                   b_habitat_high * (sq.habitat == "high_risk") +# b_habitat_half * (sq.habitat == "half_risk") +
                   sq.cost))
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings,functionality)
  P = apollo_panelProd(P,apollo_inputs,functionality)
  P = apollo_avgInterDraws(P,apollo_inputs,functionality)
  P = apollo_prepareProb(P,apollo_inputs,functionality)
  return(P)
}

#estimate model
mxl_wtp <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show output
apollo_modelOutput(mxl_wtp,
                   modelOutput_settings=list(printPVal=1))
#write output
apollo_saveOutput(mxl_wtp,
                  saveOutput_settings=list(printPVal=1))

######################
# 6. Mixed logit with all attribute interactions in WTP space
######################

apollo_control <- list(
  modelName="glaciers_CE_mxl_wtp_inter",
  modelDescr="Discrete choice experiment into preferences for glacier ES in China",
  indivID="ID",
  mixing = TRUE,
  nCores = 7,
  seed = 2104
)

# define coefficients to be estimated and starting values for each
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = 0,
                 # attribute level means
                 b_climate_1.8 = 0,
                 m_b_climate_1.5 = 0.01,
                 m_b_climate_1.2 = 0.01,
                 b_flood_freq = 0,
                 m_b_flood_10 = 0.01,
                 m_b_flood_20 = 0.01,
                 b_drought_freq = 0,
                 m_b_drought_10 = 0.01,
                 m_b_drought_20 = 0.01,
                 b_water_30 = 0,
                 m_b_water_35 = 0.01,
                 m_b_water_40 = 0.01,
                 b_habitat_high = 0,
                 m_b_habitat_half = 0.01,
                 m_log_b_cost = -1,
                 # attribute level standard deviations
                 sigma_b_climate_1.5 = 0.01,
                 sigma_b_climate_1.2 = 0.01,
                 sigma_b_flood_10 = 0.01,
                 sigma_b_flood_20 = 0.01,
                 sigma_b_drought_10 = 0.01,
                 sigma_b_drought_20 = 0.01,
                 sigma_b_water_35 = 0.01,
                 sigma_b_water_40 = 0.01,
                 sigma_b_habitat_half = 0.01,
                 sigma_log_b_cost = 1,
                 # regional interactions
                 b_basinds_climate15 = 0,
                 b_basinds_climate12 = 0,
                 b_basinds_flood10 = 0,
                 b_basinds_flood20 = 0,
                 b_basinds_drought10 = 0,
                 b_basinds_drought20 = 0,
                 b_basinds_water35 = 0,
                 b_basinds_water40 = 0,
                 b_basinds_habitathalf = 0,
                 b_basinwt_climate15 = 0,
                 b_basinwt_climate12 = 0,
                 b_basinwt_flood10 = 0,
                 b_basinwt_flood20 = 0,
                 b_basinwt_drought10 = 0,
                 b_basinwt_drought20 = 0,
                 b_basinwt_water35 = 0,
                 b_basinwt_water40 = 0,
                 b_basinwt_habitathalf = 0,
                 # socio-demographic interactions
                 b_male_climate15 = 0,
                 b_male_climate12 = 0,
                 b_male_flood10 = 0,
                 b_male_flood20 = 0,
                 b_male_drought10 = 0,
                 b_male_drought20 = 0,
                 b_male_water35 = 0,
                 b_male_water40 = 0,
                 b_male_habitathalf = 0,
                 b_urban_climate15 = 0,
                 b_urban_climate12 = 0,
                 b_urban_flood10 = 0,
                 b_urban_flood20 = 0,
                 b_urban_drought10 = 0,
                 b_urban_drought20 = 0,
                 b_urban_water35 = 0,
                 b_urban_water40 = 0,
                 b_urban_habitathalf = 0,
                 b_edu_climate15 = 0,
                 b_edu_climate12 = 0,
                 b_edu_flood10 = 0,
                 b_edu_flood20 = 0,
                 b_edu_drought10 = 0,
                 b_edu_drought20 = 0,
                 b_edu_water35 = 0,
                 b_edu_water40 = 0,
                 b_edu_habitathalf = 0,
                 b_income_climate15 = 0,
                 b_income_climate12 = 0,
                 b_income_flood10 = 0,
                 b_income_flood20 = 0,
                 b_income_drought10 = 0,
                 b_income_drought20 = 0,
                 b_income_water35 = 0,
                 b_income_water40 = 0,
                 b_income_habitathalf = 0,
                 # experience interactions
                 b_exp_drought10 = 0,
                 b_exp_drought20 = 0,
                 b_exp_flood10 = 0,
                 b_exp_flood20 = 0
                 )

# set fixed 0 coefficients
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high"
)

#define parameters for the simulation
apollo_draws = list(interDrawsType = "sobol",
                    interNDraws = 1000,
                    interNormDraws = c("draws_cost_inter",
                                       "draws_climate15_inter",
                                       "draws_climate12_inter",
                                       "draws_flood10_inter",
                                       "draws_flood20_inter",
                                       "draws_drought10_inter",
                                       "draws_drought20_inter",
                                       "draws_water35_inter",
                                       "draws_water40_inter",
                                       "draws_habitatha_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_climate_1.5"]] = m_b_climate_1.5 + sigma_b_climate_1.5 * draws_climate15_inter
  randcoeff[["b_climate_1.2"]] = m_b_climate_1.2 + sigma_b_climate_1.2 * draws_climate12_inter
  randcoeff[["b_flood_10"]] = m_b_flood_10 + sigma_b_flood_10 * draws_flood10_inter
  randcoeff[["b_flood_20"]] = m_b_flood_20 + sigma_b_flood_20 * draws_flood20_inter
  randcoeff[["b_drought_10"]] = m_b_drought_10 + sigma_b_drought_10 * draws_drought10_inter
  randcoeff[["b_drought_20"]] = m_b_drought_20 + sigma_b_drought_20 * draws_drought20_inter
  randcoeff[["b_water_35"]] = m_b_water_35 + sigma_b_water_35 * draws_water35_inter
  randcoeff[["b_water_40"]] = m_b_water_40 + sigma_b_water_40 * draws_water40_inter
  randcoeff[["b_habitat_half"]] = m_b_habitat_half + sigma_b_habitat_half * draws_habitatha_inter
  randcoeff[["b_cost"]] = -exp(m_log_b_cost + sigma_log_b_cost * draws_cost_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  
  # define interaction terms
  b_climate15 = b_climate_1.5 + b_basinds_climate15 * (basin == "ds") + b_basinwt_climate15 * (basin == "wt") +
    b_male_climate15 * (gender == 1) + b_urban_climate15 * (household == 1) + b_edu_climate15 * edu + b_income_climate15 * income
  b_climate12 = b_climate_1.2 + b_basinds_climate12 * (basin == "ds") + b_basinwt_climate12 * (basin == "wt") +
    b_male_climate12 * (gender == 1) + b_urban_climate12 * (household == 1) + b_edu_climate12 * edu + b_income_climate12 * income
  b_flood10 = b_flood_10 + b_basinds_flood10 * (basin == "ds") + b_basinwt_flood10 * (basin == "wt") +
    b_male_flood10 * (gender == 1) + b_urban_flood10 * (household == 1) + b_edu_flood10 * edu + b_income_flood10 * income +
    b_exp_flood10 * exp_flood
  b_flood20 = b_flood_20 + b_basinds_flood20 * (basin == "ds") + b_basinwt_flood20 * (basin == "wt") +
    b_male_flood20 * (gender == 1) + b_urban_flood20 * (household == 1) + b_edu_flood20 * edu + b_income_flood20 * income +
    b_exp_flood20 * exp_flood
  b_drought10 = b_drought_10 + b_basinds_drought10 * (basin == "ds") + b_basinwt_drought10 * (basin == "wt") +
    b_male_drought10 * (gender == 1) + b_urban_drought10 * (household == 1) + b_edu_drought10 * edu + b_income_drought10 * income +
    b_exp_drought10 * exp_drought
  b_drought20 = b_drought_20 + b_basinds_drought20 * (basin == "ds") + b_basinwt_drought20 * (basin == "wt") +
    b_male_drought20 * (gender == 1) + b_urban_drought20 * (household == 1) + b_edu_drought20 * edu + b_income_drought20 * income +
    b_exp_drought20 * exp_drought
  b_water35 = b_water_35 + b_basinds_water35 * (basin == "ds") + b_basinwt_water35 * (basin == "wt") +
    b_male_water35 * (gender == 1) + b_urban_water35 * (household == 1) + b_edu_water35 * edu + b_income_water35 * income
  b_water40 = b_water_40 + b_basinds_water40 * (basin == "ds") + b_basinwt_water40 * (basin == "wt") +
    b_male_water40 * (gender == 1) + b_urban_water40 * (household == 1) + b_edu_water40 * edu + b_income_water40 * income
  b_habitathalf = b_habitat_half + b_basinds_habitathalf * (basin == "ds") + b_basinwt_habitathalf * (basin == "wt") +
    b_male_habitathalf * (gender == 1) + b_urban_habitathalf * (household == 1) + b_edu_habitathalf * edu + b_income_habitathalf * income
  
  V = list()
  
  V[['alt1']] = (asc_1 + b_cost * (
    b_climate_1.8 * (alt1.climate == "1.8_degrees") + b_climate15 * (alt1.climate == "1.5_degrees") + b_climate12 * (alt1.climate == "1.2_degrees") +
      b_flood_freq * (alt1.flood == "frequent") + b_flood10 * (alt1.flood == "10_%") + b_flood20 * (alt1.flood == "20_%") +
      b_drought_freq * (alt1.drought == "frequent") + b_drought10 * (alt1.drought == "10_%") + b_drought20 * (alt1.drought == "20_%") + 
      b_water_30 * (alt1.freshwater == "30_mil") + b_water35 * (alt1.freshwater == "35_mil") + b_water40 * (alt1.freshwater == "40_mil") +
      b_habitat_high * (alt1.habitat == "high_risk") + b_habitathalf * (alt1.habitat == "half_risk") +
      alt1.cost))
  V[['alt2']] = (asc_2 + b_cost * (
    b_climate_1.8 * (alt2.climate == "1.8_degrees") + b_climate15 * (alt2.climate == "1.5_degrees") + b_climate12 * (alt2.climate == "1.2_degrees") +
      b_flood_freq * (alt2.flood == "frequent") + b_flood10 * (alt2.flood == "10_%") + b_flood20 * (alt2.flood == "20_%") +
      b_drought_freq * (alt2.drought == "frequent") + b_drought10 * (alt2.drought == "10_%") + b_drought20 * (alt2.drought == "20_%") + 
      b_water_30 * (alt2.freshwater == "30_mil") + b_water35 * (alt2.freshwater == "35_mil") + b_water40 * (alt2.freshwater == "40_mil") +
      b_habitat_high * (alt2.habitat == "high_risk") + b_habitathalf * (alt2.habitat == "half_risk") +
      alt2.cost))
  V[['alt3']] = (asc_3 + b_cost* (
    b_climate_1.8 * (sq.climate == "1.8_degrees") +
      b_flood_freq * (sq.flood == "frequent") +
      b_drought_freq * (sq.drought == "frequent") +
      b_water_30 * (sq.freshwater == "30_mil") +
      b_habitat_high * (sq.habitat == "high_risk") +
      sq.cost))
  
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

#estimate model
mxl_wtp_inter <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show outputs
apollo_modelOutput(mxl_wtp_inter,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mxl_wtp_inter,
                  saveOutput_settings=list(printPVal=1))

######################
# 7. Mixed logit with SQ interactions in WTP space
######################

apollo_control <- list(
  modelName="glaciers_CE_mxl_wtp_inter_sq",
  modelDescr="Discrete choice experiment into preferences for glacier ES in China",
  indivID="ID",
  mixing = TRUE,
  nCores = 7,
  seed = 2104
)

# define coefficients to be estimated and starting values for each
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = 0,
                 # attribute level means
                 b_climate_1.8 = 0,
                 m_b_climate_1.5 = 0.01,
                 m_b_climate_1.2 = 0.01,
                 b_flood_freq = 0,
                 m_b_flood_10 = 0.01,
                 m_b_flood_20 = 0.01,
                 b_drought_freq = 0,
                 m_b_drought_10 = 0.01,
                 m_b_drought_20 = 0.01,
                 b_water_30 = 0,
                 m_b_water_35 = 0.01,
                 m_b_water_40 = 0.01,
                 b_habitat_high = 0,
                 m_b_habitat_half = 0.01,
                 m_log_b_cost = -1,
                 # attribute level standard deviations
                 sigma_b_climate_1.5 = 0.01,
                 sigma_b_climate_1.2 = 0.01,
                 sigma_b_flood_10 = 0.01,
                 sigma_b_flood_20 = 0.01,
                 sigma_b_drought_10 = 0.01,
                 sigma_b_drought_20 = 0.01,
                 sigma_b_water_35 = 0.01,
                 sigma_b_water_40 = 0.01,
                 sigma_b_habitat_half = 0.01,
                 sigma_log_b_cost = 1,
                 # SQ interactions
                 sq_disappear = 0,
                 #sq_way_climate = 0,
                 sq_way_pollution = 0,
                 sq_way_park = 0
)

# set fixed 0 coefficients
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high"
)

#define parameters for the simulation
apollo_draws = list(interDrawsType = "sobol",
                    interNDraws = 1000,
                    interNormDraws = c("draws_cost_inter",
                                       "draws_climate15_inter",
                                       "draws_climate12_inter",
                                       "draws_flood10_inter",
                                       "draws_flood20_inter",
                                       "draws_drought10_inter",
                                       "draws_drought20_inter",
                                       "draws_water35_inter",
                                       "draws_water40_inter",
                                       "draws_habitatha_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_climate_1.5"]] = m_b_climate_1.5 + sigma_b_climate_1.5 * draws_climate15_inter
  randcoeff[["b_climate_1.2"]] = m_b_climate_1.2 + sigma_b_climate_1.2 * draws_climate12_inter
  randcoeff[["b_flood_10"]] = m_b_flood_10 + sigma_b_flood_10 * draws_flood10_inter
  randcoeff[["b_flood_20"]] = m_b_flood_20 + sigma_b_flood_20 * draws_flood20_inter
  randcoeff[["b_drought_10"]] = m_b_drought_10 + sigma_b_drought_10 * draws_drought10_inter
  randcoeff[["b_drought_20"]] = m_b_drought_20 + sigma_b_drought_20 * draws_drought20_inter
  randcoeff[["b_water_35"]] = m_b_water_35 + sigma_b_water_35 * draws_water35_inter
  randcoeff[["b_water_40"]] = m_b_water_40 + sigma_b_water_40 * draws_water40_inter
  randcoeff[["b_habitat_half"]] = m_b_habitat_half + sigma_b_habitat_half * draws_habitatha_inter
  randcoeff[["b_cost"]] = -exp(m_log_b_cost + sigma_log_b_cost * draws_cost_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  
  # define interaction terms
  asc_3 + sq_disappear * disappear + sq_way_pollution * way_pollution.control + sq_way_park * way_national.park #+ sq_way_climate
  
  V = list()
  
  V[['alt1']] = (asc_1 + b_cost * (
    b_climate_1.8 * (alt1.climate == "1.8_degrees") + b_climate_1.5 * (alt1.climate == "1.5_degrees") + b_climate_1.2 * (alt1.climate == "1.2_degrees") +
      b_flood_freq * (alt1.flood == "frequent") + b_flood_10 * (alt1.flood == "10_%") + b_flood_20 * (alt1.flood == "20_%") +
      b_drought_freq * (alt1.drought == "frequent") + b_drought_10 * (alt1.drought == "10_%") + b_drought_20 * (alt1.drought == "20_%") + 
      b_water_30 * (alt1.freshwater == "30_mil") + b_water_35 * (alt1.freshwater == "35_mil") + b_water_40 * (alt1.freshwater == "40_mil") +
      b_habitat_high * (alt1.habitat == "high_risk") + b_habitat_half * (alt1.habitat == "half_risk") +
      alt1.cost))
  V[['alt2']] = (asc_2 + b_cost * (
    b_climate_1.8 * (alt2.climate == "1.8_degrees") + b_climate_1.5 * (alt2.climate == "1.5_degrees") + b_climate_1.2 * (alt2.climate == "1.2_degrees") +
      b_flood_freq * (alt2.flood == "frequent") + b_flood_10 * (alt2.flood == "10_%") + b_flood_20 * (alt2.flood == "20_%") +
      b_drought_freq * (alt2.drought == "frequent") + b_drought_10 * (alt2.drought == "10_%") + b_drought_20 * (alt2.drought == "20_%") + 
      b_water_30 * (alt2.freshwater == "30_mil") + b_water_35 * (alt2.freshwater == "35_mil") + b_water_40 * (alt2.freshwater == "40_mil") +
      b_habitat_high * (alt2.habitat == "high_risk") + b_habitat_half * (alt2.habitat == "half_risk") +
      alt2.cost))
  V[['alt3']] = (asc_sq + b_cost* (
    b_climate_1.8 * (sq.climate == "1.8_degrees") +
      b_flood_freq * (sq.flood == "frequent") +
      b_drought_freq * (sq.drought == "frequent") +
      b_water_30 * (sq.freshwater == "30_mil") +
      b_habitat_high * (sq.habitat == "high_risk") +
      sq.cost))
  
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

#estimate model
mxl_wtp_sq <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show outputs
apollo_modelOutput(mxl_wtp_sq,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mxl_wtp_sq,
                  saveOutput_settings=list(printPVal=1))

######################
# 8. Mixed logit with all relevant interactions in preference space
######################

apollo_control <- list(
  modelName = "glaciers_CE_mxl_inter",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  mixing = TRUE,
  nCores = 20,
  seed = 2104
)

# define coefficients to be estimated and starting values for each
# starting values based on simple mixed logit and some intuitions
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = -1,
                 # attribute level means
                 b_climate_1.8 = 0,
                 m_b_climate_1.5 = 0.1,
                 m_b_climate_1.2 = 0.1,
                 b_flood_freq = 0,
                 m_b_flood_10 = 0.1,
                 m_b_flood_20 = 0.1,
                 b_drought_freq = 0,
                 m_b_drought_10 = 0.1,
                 m_b_drought_20 = 0.1,
                 b_water_30 = 0,
                 m_b_water_35 = 0,
                 m_b_water_40 = 0.01,
                 b_habitat_high = 0,
                 m_b_habitat_half = 0.1,
                 m_log_b_cost = -5,
                 # attribute level standard deviations
                 sigma_b_climate_1.5 = 0.01,
                 sigma_b_climate_1.2 = 0.01,
                 sigma_b_flood_10 = 0.01,
                 sigma_b_flood_20 = 0.01,
                 sigma_b_drought_10 = 0.01,
                 sigma_b_drought_20 = 0.01,
                 sigma_b_water_35 = 0.01,
                 sigma_b_water_40 = 0.01,
                 sigma_b_habitat_half = 0.01,
                 sigma_log_b_cost = 1,
                 # regional interactions
                 b_basinds_climate15 = 0,
                 b_basinds_climate12 = 0,
                 b_basinds_flood10 = 0,
                 b_basinds_flood20 = 0,
                 b_basinds_drought10 = 0,
                 b_basinds_drought20 = 0,
                 b_basinds_water35 = 0,
                 b_basinds_water40 = 0,
                 b_basinds_habitathalf = 0,
                 b_basinwt_climate15 = 0,
                 b_basinwt_climate12 = 0,
                 b_basinwt_flood10 = 0,
                 b_basinwt_flood20 = 0,
                 b_basinwt_drought10 = 0,
                 b_basinwt_drought20 = 0,
                 b_basinwt_water35 = 0,
                 b_basinwt_water40 = 0,
                 b_basinwt_habitathalf = 0,
                 # socio-demographic interactions
                 b_male_climate15 = 0,
                 b_male_climate12 = 0,
                 b_male_flood10 = 0,
                 b_male_flood20 = 0,
                 b_male_drought10 = 0,
                 b_male_drought20 = 0,
                 b_male_water35 = 0,
                 b_male_water40 = 0,
                 b_male_habitathalf = 0,
                 b_urban_climate15 = 0,
                 b_urban_climate12 = 0,
                 b_urban_flood10 = 0,
                 b_urban_flood20 = 0,
                 b_urban_drought10 = 0,
                 b_urban_drought20 = 0,
                 b_urban_water35 = 0,
                 b_urban_water40 = 0,
                 b_urban_habitathalf = 0,
                 b_edu_climate15 = 0,
                 b_edu_climate12 = 0,
                 b_edu_flood10 = 0,
                 b_edu_flood20 = 0,
                 b_edu_drought10 = 0,
                 b_edu_drought20 = 0,
                 b_edu_water35 = 0,
                 b_edu_water40 = 0,
                 b_edu_habitathalf = 0,
                 b_income_climate15 = 0,
                 b_income_climate12 = 0,
                 b_income_flood10 = 0,
                 b_income_flood20 = 0,
                 b_income_drought10 = 0,
                 b_income_drought20 = 0,
                 b_income_water35 = 0,
                 b_income_water40 = 0,
                 b_income_habitathalf = 0,
                 # experience interactions
                 b_exp_drought10 = 0,
                 b_exp_drought20 = 0,
                 b_exp_flood10 = 0,
                 b_exp_flood20 = 0,
                 # SQ interactions
                 sq_disappear = 0.001,
                 sq_way_pollution = 0,
                 sq_way_park = 0
)

# set fixed 0 coefficients
apollo_fixed <- c("asc_2",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high"
)

#define parameters for the simulation
apollo_draws = list(interDrawsType = "sobol",
                    interNDraws = 1000,
                    interNormDraws = c("draws_cost_inter",
                                       "draws_climate15_inter",
                                       "draws_climate12_inter",
                                       "draws_flood10_inter",
                                       "draws_flood20_inter",
                                       "draws_drought10_inter",
                                       "draws_drought20_inter",
                                       "draws_water35_inter",
                                       "draws_water40_inter",
                                       "draws_habitatha_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_climate_1.5"]] = m_b_climate_1.5 + sigma_b_climate_1.5 * draws_climate15_inter
  randcoeff[["b_climate_1.2"]] = m_b_climate_1.2 + sigma_b_climate_1.2 * draws_climate12_inter
  randcoeff[["b_flood_10"]] = m_b_flood_10 + sigma_b_flood_10 * draws_flood10_inter
  randcoeff[["b_flood_20"]] = m_b_flood_20 + sigma_b_flood_20 * draws_flood20_inter
  randcoeff[["b_drought_10"]] = m_b_drought_10 + sigma_b_drought_10 * draws_drought10_inter
  randcoeff[["b_drought_20"]] = m_b_drought_20 + sigma_b_drought_20 * draws_drought20_inter
  randcoeff[["b_water_35"]] = m_b_water_35 + sigma_b_water_35 * draws_water35_inter
  randcoeff[["b_water_40"]] = m_b_water_40 + sigma_b_water_40 * draws_water40_inter
  randcoeff[["b_habitat_half"]] = m_b_habitat_half + sigma_b_habitat_half * draws_habitatha_inter
  randcoeff[["b_cost"]] = -exp(m_log_b_cost + sigma_log_b_cost * draws_cost_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  
  # define interaction terms (income log'ed to avoid scale issues in estimation)
  asc_sq = asc_3 + sq_disappear * disappear + sq_way_pollution * way_pollution.control + sq_way_park * way_national.park
  b_climate15 = b_climate_1.5 + b_basinds_climate15 * (basin == "ds") + b_basinwt_climate15 * (basin == "wt") +
    b_male_climate15 * (gender == 1) + b_urban_climate15 * (household == 1) + b_edu_climate15 * edu + b_income_climate15 * log(income)
  b_climate12 = b_climate_1.2 + b_basinds_climate12 * (basin == "ds") + b_basinwt_climate12 * (basin == "wt") +
    b_male_climate12 * (gender == 1) + b_urban_climate12 * (household == 1) + b_edu_climate12 * edu + b_income_climate12 * log(income)
  b_flood10 = b_flood_10 + b_basinds_flood10 * (basin == "ds") + b_basinwt_flood10 * (basin == "wt") +
    b_male_flood10 * (gender == 1) + b_urban_flood10 * (household == 1) + b_edu_flood10 * edu + b_income_flood10 * log(income) +
    b_exp_flood10 * exp_flood
  b_flood20 = b_flood_20 + b_basinds_flood20 * (basin == "ds") + b_basinwt_flood20 * (basin == "wt") +
    b_male_flood20 * (gender == 1) + b_urban_flood20 * (household == 1) + b_edu_flood20 * edu + b_income_flood20 * log(income) +
    b_exp_flood20 * exp_flood
  b_drought10 = b_drought_10 + b_basinds_drought10 * (basin == "ds") + b_basinwt_drought10 * (basin == "wt") +
    b_male_drought10 * (gender == 1) + b_urban_drought10 * (household == 1) + b_edu_drought10 * edu + b_income_drought10 * log(income) +
    b_exp_drought10 * exp_drought
  b_drought20 = b_drought_20 + b_basinds_drought20 * (basin == "ds") + b_basinwt_drought20 * (basin == "wt") +
    b_male_drought20 * (gender == 1) + b_urban_drought20 * (household == 1) + b_edu_drought20 * edu + b_income_drought20 * log(income) +
    b_exp_drought20 * exp_drought
  b_water35 = b_water_35 + b_basinds_water35 * (basin == "ds") + b_basinwt_water35 * (basin == "wt") +
    b_male_water35 * (gender == 1) + b_urban_water35 * (household == 1) + b_edu_water35 * edu + b_income_water35 * log(income)
  b_water40 = b_water_40 + b_basinds_water40 * (basin == "ds") + b_basinwt_water40 * (basin == "wt") +
    b_male_water40 * (gender == 1) + b_urban_water40 * (household == 1) + b_edu_water40 * edu + b_income_water40 * log(income)
  b_habitathalf = b_habitat_half + b_basinds_habitathalf * (basin == "ds") + b_basinwt_habitathalf * (basin == "wt") +
    b_male_habitathalf * (gender == 1) + b_urban_habitathalf * (household == 1) + b_edu_habitathalf * edu + b_income_habitathalf * log(income)
  
  V = list()
  
  V[['alt1']] = (asc_1 + b_climate_1.8 * (alt1.climate == "1.8_degrees") + b_climate15 * (alt1.climate == "1.5_degrees") + b_climate12 * (alt1.climate == "1.2_degrees") +
      b_flood_freq * (alt1.flood == "frequent") + b_flood10 * (alt1.flood == "10_%") + b_flood20 * (alt1.flood == "20_%") +
      b_drought_freq * (alt1.drought == "frequent") + b_drought10 * (alt1.drought == "10_%") + b_drought20 * (alt1.drought == "20_%") + 
      b_water_30 * (alt1.freshwater == "30_mil") + b_water35 * (alt1.freshwater == "35_mil") + b_water40 * (alt1.freshwater == "40_mil") +
      b_habitat_high * (alt1.habitat == "high_risk") + b_habitathalf * (alt1.habitat == "half_risk") +
      b_cost * alt1.cost)
  V[['alt2']] = (asc_2 + b_climate_1.8 * (alt2.climate == "1.8_degrees") + b_climate15 * (alt2.climate == "1.5_degrees") + b_climate12 * (alt2.climate == "1.2_degrees") +
      b_flood_freq * (alt2.flood == "frequent") + b_flood10 * (alt2.flood == "10_%") + b_flood20 * (alt2.flood == "20_%") +
      b_drought_freq * (alt2.drought == "frequent") + b_drought10 * (alt2.drought == "10_%") + b_drought20 * (alt2.drought == "20_%") + 
      b_water_30 * (alt2.freshwater == "30_mil") + b_water35 * (alt2.freshwater == "35_mil") + b_water40 * (alt2.freshwater == "40_mil") +
      b_habitat_high * (alt2.habitat == "high_risk") + b_habitathalf * (alt2.habitat == "half_risk") +
      b_cost * alt2.cost)
  V[['alt3']] = (asc_sq + b_climate_1.8 * (sq.climate == "1.8_degrees") +
      b_flood_freq * (sq.flood == "frequent") +
      b_drought_freq * (sq.drought == "frequent") +
      b_water_30 * (sq.freshwater == "30_mil") +
      b_habitat_high * (sq.habitat == "high_risk") +
      b_cost * sq.cost)
  
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_avgInterDraws(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

#estimate model
mxl_all_inter <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show outputs
apollo_modelOutput(mxl_all_inter,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mxl_all_inter,
                  saveOutput_settings=list(printPVal=1))

######################
# 9. Multinomial logit with all relevant interactions in preference space
######################

apollo_control <- list(
  modelName = "glaciers_CE_mnl_inter",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID"
)

# define coefficients to be estimated and starting values for each
# starting values based on simple mixed logit and some intuitions
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = 0,
                 # attribute levels
                 b_climate_1.8 = 0,
                 b_climate_1.5 = 0.1,
                 b_climate_1.2 = 0.1,
                 b_flood_freq = 0,
                 b_flood_10 = 0.1,
                 b_flood_20 = 0.1,
                 b_drought_freq = 0,
                 b_drought_10 = 0.1,
                 b_drought_20 = 0.1,
                 b_water_30 = 0,
                 b_water_35 = 0,
                 b_water_40 = 0.01,
                 b_habitat_high = 0,
                 b_habitat_half = 0.1,
                 b_cost = -1,
                 # regional interactions
                 b_basinds_climate15 = 0,
                 b_basinds_climate12 = 0,
                 b_basinds_flood10 = 0,
                 b_basinds_flood20 = 0,
                 b_basinds_drought10 = 0,
                 b_basinds_drought20 = 0,
                 b_basinds_water35 = 0,
                 b_basinds_water40 = 0,
                 b_basinds_habitathalf = 0,
                 b_basinwt_climate15 = 0,
                 b_basinwt_climate12 = 0,
                 b_basinwt_flood10 = 0,
                 b_basinwt_flood20 = 0,
                 b_basinwt_drought10 = 0,
                 b_basinwt_drought20 = 0,
                 b_basinwt_water35 = 0,
                 b_basinwt_water40 = 0,
                 b_basinwt_habitathalf = 0,
                 # socio-demographic interactions
                 b_male_climate15 = 0,
                 b_male_climate12 = 0,
                 b_male_flood10 = 0,
                 b_male_flood20 = 0,
                 b_male_drought10 = 0,
                 b_male_drought20 = 0,
                 b_male_water35 = 0,
                 b_male_water40 = 0,
                 b_male_habitathalf = 0,
                 b_urban_climate15 = 0,
                 b_urban_climate12 = 0,
                 b_urban_flood10 = 0,
                 b_urban_flood20 = 0,
                 b_urban_drought10 = 0,
                 b_urban_drought20 = 0,
                 b_urban_water35 = 0,
                 b_urban_water40 = 0,
                 b_urban_habitathalf = 0,
                 b_edu_climate15 = 0,
                 b_edu_climate12 = 0,
                 b_edu_flood10 = 0,
                 b_edu_flood20 = 0,
                 b_edu_drought10 = 0,
                 b_edu_drought20 = 0,
                 b_edu_water35 = 0,
                 b_edu_water40 = 0,
                 b_edu_habitathalf = 0,
                 b_income_climate15 = 0,
                 b_income_climate12 = 0,
                 b_income_flood10 = 0,
                 b_income_flood20 = 0,
                 b_income_drought10 = 0,
                 b_income_drought20 = 0,
                 b_income_water35 = 0,
                 b_income_water40 = 0,
                 b_income_habitathalf = 0,
                 # experience interactions
                 b_exp_drought10 = 0,
                 b_exp_drought20 = 0,
                 b_exp_flood10 = 0,
                 b_exp_flood20 = 0,
                 # SQ interactions
                 sq_disappear = 0.001,
                 sq_way_pollution = 0,
                 sq_way_park = 0
)

# set fixed 0 coefficients
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high"
)

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  
  # define interaction terms
  asc_sq = asc_3 + sq_disappear * disappear + sq_way_pollution * way_pollution.control + sq_way_park * way_national.park
  b_climate15 = b_climate_1.5 + b_basinds_climate15 * (basin == "ds") + b_basinwt_climate15 * (basin == "wt") +
    b_male_climate15 * (gender == 1) + b_urban_climate15 * (household == 1) + b_edu_climate15 * edu + b_income_climate15 * income
  b_climate12 = b_climate_1.2 + b_basinds_climate12 * (basin == "ds") + b_basinwt_climate12 * (basin == "wt") +
    b_male_climate12 * (gender == 1) + b_urban_climate12 * (household == 1) + b_edu_climate12 * edu + b_income_climate12 * income
  b_flood10 = b_flood_10 + b_basinds_flood10 * (basin == "ds") + b_basinwt_flood10 * (basin == "wt") +
    b_male_flood10 * (gender == 1) + b_urban_flood10 * (household == 1) + b_edu_flood10 * edu + b_income_flood10 * income +
    b_exp_flood10 * exp_flood
  b_flood20 = b_flood_20 + b_basinds_flood20 * (basin == "ds") + b_basinwt_flood20 * (basin == "wt") +
    b_male_flood20 * (gender == 1) + b_urban_flood20 * (household == 1) + b_edu_flood20 * edu + b_income_flood20 * income +
    b_exp_flood20 * exp_flood
  b_drought10 = b_drought_10 + b_basinds_drought10 * (basin == "ds") + b_basinwt_drought10 * (basin == "wt") +
    b_male_drought10 * (gender == 1) + b_urban_drought10 * (household == 1) + b_edu_drought10 * edu + b_income_drought10 * income +
    b_exp_drought10 * exp_drought
  b_drought20 = b_drought_20 + b_basinds_drought20 * (basin == "ds") + b_basinwt_drought20 * (basin == "wt") +
    b_male_drought20 * (gender == 1) + b_urban_drought20 * (household == 1) + b_edu_drought20 * edu + b_income_drought20 * income +
    b_exp_drought20 * exp_drought
  b_water35 = b_water_35 + b_basinds_water35 * (basin == "ds") + b_basinwt_water35 * (basin == "wt") +
    b_male_water35 * (gender == 1) + b_urban_water35 * (household == 1) + b_edu_water35 * edu + b_income_water35 * income
  b_water40 = b_water_40 + b_basinds_water40 * (basin == "ds") + b_basinwt_water40 * (basin == "wt") +
    b_male_water40 * (gender == 1) + b_urban_water40 * (household == 1) + b_edu_water40 * edu + b_income_water40 * income
  b_habitathalf = b_habitat_half + b_basinds_habitathalf * (basin == "ds") + b_basinwt_habitathalf * (basin == "wt") +
    b_male_habitathalf * (gender == 1) + b_urban_habitathalf * (household == 1) + b_edu_habitathalf * edu + b_income_habitathalf * income
  
  V = list()
  
  V[['alt1']] = (asc_1 + b_climate_1.8 * (alt1.climate == "1.8_degrees") + b_climate15 * (alt1.climate == "1.5_degrees") + b_climate12 * (alt1.climate == "1.2_degrees") +
                   b_flood_freq * (alt1.flood == "frequent") + b_flood10 * (alt1.flood == "10_%") + b_flood20 * (alt1.flood == "20_%") +
                   b_drought_freq * (alt1.drought == "frequent") + b_drought10 * (alt1.drought == "10_%") + b_drought20 * (alt1.drought == "20_%") + 
                   b_water_30 * (alt1.freshwater == "30_mil") + b_water35 * (alt1.freshwater == "35_mil") + b_water40 * (alt1.freshwater == "40_mil") +
                   b_habitat_high * (alt1.habitat == "high_risk") + b_habitathalf * (alt1.habitat == "half_risk") +
                   b_cost * alt1.cost)
  V[['alt2']] = (asc_2 + b_climate_1.8 * (alt2.climate == "1.8_degrees") + b_climate15 * (alt2.climate == "1.5_degrees") + b_climate12 * (alt2.climate == "1.2_degrees") +
                   b_flood_freq * (alt2.flood == "frequent") + b_flood10 * (alt2.flood == "10_%") + b_flood20 * (alt2.flood == "20_%") +
                   b_drought_freq * (alt2.drought == "frequent") + b_drought10 * (alt2.drought == "10_%") + b_drought20 * (alt2.drought == "20_%") + 
                   b_water_30 * (alt2.freshwater == "30_mil") + b_water35 * (alt2.freshwater == "35_mil") + b_water40 * (alt2.freshwater == "40_mil") +
                   b_habitat_high * (alt2.habitat == "high_risk") + b_habitathalf * (alt2.habitat == "half_risk") +
                   b_cost * alt2.cost)
  V[['alt3']] = (asc_sq + b_climate_1.8 * (sq.climate == "1.8_degrees") +
                   b_flood_freq * (sq.flood == "frequent") +
                   b_drought_freq * (sq.drought == "frequent") +
                   b_water_30 * (sq.freshwater == "30_mil") +
                   b_habitat_high * (sq.habitat == "high_risk") +
                   b_cost * sq.cost)
  
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

#estimate model
mnl_all_inter <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show outputs
apollo_modelOutput(mnl_all_inter,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mnl_all_inter,
                  saveOutput_settings=list(printPVal=1))

######################
# 10. Mixed logit with regional interactions in WTP space
######################

apollo_control <- list(
  modelName="glaciers_CE_mxl_wtp_reg",
  modelDescr="Discrete choice experiment into preferences for glacier ES in China",
  indivID="ID",
  mixing = TRUE,
  nCores = 7,
  seed = 2104
)

# define coefficients to be estimated and starting values for each
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = 0,
                 b_climate_1.8 = 0,
                 m_b_climate_1.5 = 0.01,
                 m_b_climate_1.2 = 0.01,
                 b_flood_freq = 0,
                 m_b_flood_10 = 0.01,
                 m_b_flood_20 = 0.01,
                 b_drought_freq = 0,
                 m_b_drought_10 = 0.01,
                 m_b_drought_20 = 0.01,
                 b_water_30 = 0,
                 m_b_water_35 = 0.01,
                 m_b_water_40 = 0.01,
                 b_habitat_high = 0,
                 m_b_habitat_half = 0.01,
                 m_log_b_cost = -1,
                 sigma_b_climate_1.5 = 0.01,
                 sigma_b_climate_1.2 = 0.01,
                 sigma_b_flood_10 = 0.01,
                 sigma_b_flood_20 = 0.01,
                 sigma_b_drought_10 = 0.01,
                 sigma_b_drought_20 = 0.01,
                 sigma_b_water_35 = 0.01,
                 sigma_b_water_40 = 0.01,
                 sigma_b_habitat_half = 0.01,
                 sigma_log_b_cost = 1,
                 b_basinds_climate15 = 0,
                 b_basinds_climate12 = 0,
                 b_basinds_flood10 = 0,
                 b_basinds_flood20 = 0,
                 b_basinds_drought10 = 0,
                 b_basinds_drought20 = 0,
                 b_basinds_water35 = 0,
                 b_basinds_water40 = 0,
                 b_basinds_habitathalf = 0,
                 b_basinwt_climate15 = 0,
                 b_basinwt_climate12 = 0,
                 b_basinwt_flood10 = 0,
                 b_basinwt_flood20 = 0,
                 b_basinwt_drought10 = 0,
                 b_basinwt_drought20 = 0,
                 b_basinwt_water35 = 0,
                 b_basinwt_water40 = 0,
                 b_basinwt_habitathalf = 0)

#set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high"
)

#define parameters for the simulation
apollo_draws = list(interDrawsType = "sobol",
                    interNDraws = 1000,
                    interNormDraws = c("draws_cost_inter",
                                       #"draws_climate18_inter",
                                       "draws_climate15_inter",
                                       "draws_climate12_inter",
                                       #"draws_floodfreq_inter",
                                       "draws_flood10_inter",
                                       "draws_flood20_inter",
                                       #"draws_droughtfreq_inter",
                                       "draws_drought10_inter",
                                       "draws_drought20_inter",
                                       #"draws_water30_inter",
                                       "draws_water35_inter",
                                       "draws_water40_inter",
                                       #"draws_habitathi_inter",
                                       "draws_habitatha_inter")
)

#define random coefficients
apollo_randCoeff = function(apollo_beta,apollo_inputs){
  randcoeff = list()
  randcoeff[["b_climate_1.5"]] = m_b_climate_1.5 + sigma_b_climate_1.5 * draws_climate15_inter
  randcoeff[["b_climate_1.2"]] = m_b_climate_1.2 + sigma_b_climate_1.2 * draws_climate12_inter
  randcoeff[["b_flood_10"]] = m_b_flood_10 + sigma_b_flood_10 * draws_flood10_inter
  randcoeff[["b_flood_20"]] = m_b_flood_20 + sigma_b_flood_20 * draws_flood20_inter
  randcoeff[["b_drought_10"]] = m_b_drought_10 + sigma_b_drought_10 * draws_drought10_inter
  randcoeff[["b_drought_20"]] = m_b_drought_20 + sigma_b_drought_20 * draws_drought20_inter
  randcoeff[["b_water_35"]] = m_b_water_35 + sigma_b_water_35 * draws_water35_inter
  randcoeff[["b_water_40"]] = m_b_water_40 + sigma_b_water_40 * draws_water40_inter
  randcoeff[["b_habitat_half"]] = m_b_habitat_half + sigma_b_habitat_half * draws_habitatha_inter
  randcoeff[["b_cost"]] = -exp(m_log_b_cost + sigma_log_b_cost * draws_cost_inter)
  return(randcoeff)
}

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  
  # define interaction terms
  b_climate15 = b_climate_1.5 + b_basinds_climate15 * (basin == "ds") + b_basinwt_climate15 * (basin == "wt")
  b_climate12 = b_climate_1.2 + b_basinds_climate12 * (basin == "ds") + b_basinwt_climate12 * (basin == "wt")
  b_flood10 = b_flood_10 + b_basinds_flood10 * (basin == "ds") + b_basinwt_flood10 * (basin == "wt")
  b_flood20 = b_flood_20 + b_basinds_flood20 * (basin == "ds") + b_basinwt_flood20 * (basin == "wt")
  b_drought10 = b_drought_10 + b_basinds_drought10 * (basin == "ds") + b_basinwt_drought10 * (basin == "wt")
  b_drought20 = b_drought_20 + b_basinds_drought20 * (basin == "ds") + b_basinwt_drought20 * (basin == "wt")
  b_water35 = b_water_35 + b_basinds_water35 * (basin == "ds") + b_basinwt_water35 * (basin == "wt")
  b_water40 = b_water_40 + b_basinds_water40 * (basin == "ds") + b_basinwt_water40 * (basin == "wt")
  b_habitathalf = b_habitat_half + b_basinds_habitathalf * (basin == "ds") + b_basinwt_habitathalf * (basin == "wt")
  
  V = list()
  
  V[['alt1']] = (asc_1 + b_cost * (
    b_climate_1.8 * (alt1.climate == "1.8_degrees") + b_climate15 * (alt1.climate == "1.5_degrees") + b_climate12 * (alt1.climate == "1.2_degrees") +
      b_flood_freq * (alt1.flood == "frequent") + b_flood10 * (alt1.flood == "10_%") + b_flood20 * (alt1.flood == "20_%") +
      b_drought_freq * (alt1.drought == "frequent") + b_drought10 * (alt1.drought == "10_%") + b_drought20 * (alt1.drought == "20_%") + 
      b_water_30 * (alt1.freshwater == "30_mil") + b_water35 * (alt1.freshwater == "35_mil") + b_water40 * (alt1.freshwater == "40_mil") +
      b_habitat_high * (alt1.habitat == "high_risk") + b_habitathalf * (alt1.habitat == "half_risk") +
      alt1.cost))
  V[['alt2']] = (asc_2 + b_cost * (
    b_climate_1.8 * (alt2.climate == "1.8_degrees") + b_climate15 * (alt2.climate == "1.5_degrees") + b_climate12 * (alt2.climate == "1.2_degrees") +
      b_flood_freq * (alt2.flood == "frequent") + b_flood10 * (alt2.flood == "10_%") + b_flood20 * (alt2.flood == "20_%") +
      b_drought_freq * (alt2.drought == "frequent") + b_drought10 * (alt2.drought == "10_%") + b_drought20 * (alt2.drought == "20_%") + 
      b_water_30 * (alt2.freshwater == "30_mil") + b_water35 * (alt2.freshwater == "35_mil") + b_water40 * (alt2.freshwater == "40_mil") +
      b_habitat_high * (alt2.habitat == "high_risk") + b_habitathalf * (alt2.habitat == "half_risk") +
      alt2.cost))
  V[['alt3']] = (asc_3 + b_cost* (
    b_climate_1.8 * (sq.climate == "1.8_degrees") +# b_climate_1.5 * (sq.climate == "1.5_degrees") + b_climate_1.2 * (sq.climate == "1.2_degrees") +
      b_flood_freq * (sq.flood == "frequent") +# b_flood_10 * (sq.flood == "10_%") + b_flood_20 * (sq.flood == "20_%") +
      b_drought_freq * (sq.drought == "frequent") +# b_drought_10 * (sq.drought == "10_%") + b_drought_20 * (sq.drought == "20_%") + 
      b_water_30 * (sq.freshwater == "30_mil") +# b_water_35 * (sq.freshwater == "35_mil") + b_water_40 * (sq.freshwater == "40_mil") +
      b_habitat_high * (sq.habitat == "high_risk") +# b_habitat_half * (sq.habitat == "half_risk") +
      sq.cost))
  
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings,functionality)
  P = apollo_panelProd(P,apollo_inputs,functionality)
  P = apollo_avgInterDraws(P,apollo_inputs,functionality)
  P = apollo_prepareProb(P,apollo_inputs,functionality)
  return(P)
}

#estimate model
mxl_wtp_reg <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show outputs
apollo_modelOutput(mxl_wtp_reg,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mxl_wtp_reg,
                  saveOutput_settings=list(printPVal=1))

#######################
# 11. [test] Multinomial logit with continuous variables
#######################

# redefine variables
database$alt1.climate <- ifelse(database$alt1.climate == "1.8_degrees", 1.8,
                                ifelse(database$alt1.climate == "1.5_degrees", 1.5,
                                       ifelse(database$alt1.climate == "1.2_degrees", 1.2,
                                              database$alt1.climate)))
database$alt1.climate <- as.numeric(database$alt1.climate)
database$alt1.drought <- ifelse(database$alt1.drought == "frequent", 0,
                                ifelse(database$alt1.drought == "10_%", 10,
                                       ifelse(database$alt1.drought == "20_%", 20,
                                              database$alt1.drought)))
database$alt1.drought <- as.numeric(database$alt1.drought)
database$alt1.flood <- ifelse(database$alt1.flood == "frequent", 0,
                              ifelse(database$alt1.flood == "10_%", 10,
                                     ifelse(database$alt1.flood == "20_%", 20,
                                            database$alt1.flood)))
database$alt1.flood <- as.numeric(database$alt1.flood)
database$alt1.freshwater <- ifelse(database$alt1.freshwater == "30_mil", 30,
                                   ifelse(database$alt1.freshwater == "35_mil", 35,
                                          ifelse(database$alt1.freshwater == "40_mil", 40,
                                                 database$alt1.freshwater)))
database$alt1.freshwater <- as.numeric(database$alt1.freshwater)
database$alt1.habitat <- ifelse(database$alt1.habitat == "high_risk", 0, 50)
database$alt1.habitat <- as.numeric(database$alt1.habitat)
database$alt2.climate <- ifelse(database$alt2.climate == "1.8_degrees", 1.8,
                                ifelse(database$alt2.climate == "1.5_degrees", 1.5,
                                       ifelse(database$alt2.climate == "1.2_degrees", 1.2,
                                              database$alt2.climate)))
database$alt2.climate <- as.numeric(database$alt2.climate)
database$alt2.drought <- ifelse(database$alt2.drought == "frequent", 0,
                                ifelse(database$alt2.drought == "10_%", 10,
                                       ifelse(database$alt2.drought == "20_%", 20,
                                              database$alt2.drought)))
database$alt2.drought <- as.numeric(database$alt2.drought)
database$alt2.flood <- ifelse(database$alt2.flood == "frequent", 0,
                              ifelse(database$alt2.flood == "10_%", 10,
                                     ifelse(database$alt2.flood == "20_%", 20,
                                            database$alt2.flood)))
database$alt2.flood <- as.numeric(database$alt2.flood)
database$alt2.freshwater <- ifelse(database$alt2.freshwater == "30_mil", 30,
                                   ifelse(database$alt2.freshwater == "35_mil", 35,
                                          ifelse(database$alt2.freshwater == "40_mil", 40,
                                                 database$alt2.freshwater)))
database$alt2.freshwater <- as.numeric(database$alt2.freshwater)
database$alt2.habitat <- ifelse(database$alt2.habitat == "high_risk", 0, 50)
database$alt2.habitat <- as.numeric(database$alt2.habitat)
database$sq.climate <- 1.8
database$sq.drought <- 0
database$sq.flood <- 0
database$sq.freshwater <- 30
database$sq.habitat <- 0


# set required parameters of project
apollo_control <- list(
  modelName="glaciers_CE_mnl_cont",
  modelDescr="Discrete choice experiment into preferences for glacier ES in China",
  indivID="ID"
)

# define coefficients to be estimated and starting values for each
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = 0,
                 b_climate = 0,
                 b_flood = 0,
                 b_drought = 0,
                 b_water = 0,
                 b_habitat = 0,
                 b_cost = 0)

#set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_3")

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  V = list()
  # define utilities for each alternative (brackets required because of the breaks)
  V[['alt1']] = (asc_1 + 
                   b_climate * alt1.climate +
                   b_flood * alt1.flood +
                   b_drought * alt1.drought + 
                   b_water * alt1.freshwater +
                   b_habitat * alt1.habitat +
                   b_cost * alt1.cost)
  V[['alt2']] = (asc_2 + 
                   b_climate * alt2.climate +
                   b_flood * alt2.flood +
                   b_drought * alt2.drought + 
                   b_water * alt2.freshwater +
                   b_habitat * alt2.habitat +
                   b_cost * alt2.cost)
  V[['alt3']] = (asc_3 + 
                   b_climate * sq.climate +
                   b_flood * sq.flood +
                   b_drought * sq.drought +
                   b_water * sq.freshwater +
                   b_habitat * sq.habitat +
                   b_cost * sq.cost)
  
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  ) #model components
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

#estimate model
mnl_cont <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show outputs
apollo_modelOutput(mnl_cont,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mnl_cont,
                  saveOutput_settings=list(printPVal=1))

#######################
# 12. Save outputs as DOCX table(s)
#######################

# optionally: load required models
mxl_all_inter <- apollo_loadModel("glaciers_CE_mxl_inter")
mxl <- apollo_loadModel("glaciers_CE_mxl")
mnl <- apollo_loadModel("glaciers_CE_mnl")

# save model as texreg-readable object
mxl_all_inter2 <- quicktexregapollo(mxl_all_inter)
mxl2 <- quicktexregapollo(mxl)
# mnl2 <- quicktexregapollo(mnl)

# write output to word table
wordreg(list(mxl2, mxl_all_inter2), "glaciers_CE_output_table.docx")

        