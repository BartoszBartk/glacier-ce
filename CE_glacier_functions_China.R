############################
# Analysis of the data from a DCE on Chinese public's preferences for glacier functions/ecosystem services
# n=3000, conducted online in July 2024
# study authors: Can Zhang, Michael Beckmann, Martin Volk, Bo Su, Bartosz Bartkowski
# code author: Bartosz Bartkowski
# contact: bartosz.bartkowski@ufz.de
# This script can be found at https://github.com/BartoszBartk/glacier-ce
# Version: 25 June 2025
############################

############## Outline ###############
# 1. Packages and data
# 2. Multinomial logit
# 3. Simple mixed logit
# 4. Nested logit
# 5. Multinomial logit with interactions in preference space
# 5a. Multinomial logit with significant interactions only
# 6. Simple mixed logit in WTP space
# 7. Mixed logit with all relevant interactions in preference space
# 8. Subregion-specific multinomial logits
# 9. Latent Class Model
# 10. Save outputs as DOCX table(s)
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
  modelName = "glaciers_CE_mnl",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  outputDirectory = "Model_outputs"
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

# set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high")

# validate inputs
apollo_inputs <- apollo_validateInputs()

# define model (for details of the required function elements, see apollo documentation)
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
    b_climate_1.8 * (sq.climate == "1.8_degrees") +
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

# estimate model
mnl <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# show outputs
apollo_modelOutput(mnl,
                   modelOutput_settings=list(printPVal=1))

# write outputs
apollo_saveOutput(mnl,
                  saveOutput_settings=list(printPVal=1))

# calculate WTP mean and standard errors
deltaMethod_settings = list(expression = c(wtp_climate15 = "-b_climate_1.5/b_cost",
                                           wtp_climate12 = "-b_climate_1.2/b_cost",
                                           wtp_flood10 = "-b_flood_10/b_cost",
                                           wtp_flood20 = "-b_flood_20/b_cost",
                                           wtp_drought10 = "-b_drought_10/b_cost",
                                           wtp_drought20 = "-b_drought_20/b_cost",
                                           wtp_water35 = "-b_water_35/b_cost",
                                           wtp_water40 = "-b_water_40/b_cost",
                                           wtp_habhalf = "-b_habitat_half/b_cost"
))
apollo_deltaMethod(mnl, deltaMethod_settings)

#######################
# 3. Simple mixed logit
#######################

# implement a simple mixed logit model to allow for inter-individual preference heterogeneity
apollo_control = list(
  modelName = "glaciers_CE_mxl",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  mixing = TRUE,
  nCores = 7,
  seed = 2104,
  outputDirectory = "Model_outputs"
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
                 sigma_log_b_cost = 1)

# set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high"
)

# define parameters for the simulation
apollo_draws = list(
  interDrawsType = "sobol",
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

# define random coefficients
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

# validate inputs
apollo_inputs <- apollo_validateInputs()

# define model (for details, see apollo documentation)
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
                   b_climate_1.8 * (sq.climate == "1.8_degrees") +
                   b_flood_freq * (sq.flood == "frequent") +
                   b_drought_freq * (sq.drought == "frequent") +
                   b_water_30 * (sq.freshwater == "30_mil") +
                   b_habitat_high * (sq.habitat == "high_risk") +
                   b_cost * sq.cost)
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  )
  P[["model"]] = apollo_mnl(mnl_settings,functionality)
  P = apollo_panelProd(P,apollo_inputs,functionality)
  P = apollo_avgInterDraws(P,apollo_inputs,functionality)
  P = apollo_prepareProb(P,apollo_inputs,functionality)
  return(P)
}

# estimate model
mxl <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# show output
apollo_modelOutput(mxl,
                   modelOutput_settings=list(printPVal=1))
# write output
apollo_saveOutput(mxl,
                  saveOutput_settings=list(printPVal=1))

# calculate WTP with help of Delta method
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
# 4. Nested logit
######################
# set required parameters of project
apollo_control <- list(
  modelName = "glaciers_CE_nl",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  outputDirectory = "Model_outputs"
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
                 b_cost = 0,
                 lambdanSQ = 1)

# set fixed 0 coefficients (here: only one of the alternative-specific constants)
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
                   b_climate_1.8 * (sq.climate == "1.8_degrees") +
                   b_flood_freq * (sq.flood == "frequent") +
                   b_drought_freq * (sq.drought == "frequent") +
                   b_water_30 * (sq.freshwater == "30_mil") +
                   b_habitat_high * (sq.habitat == "high_risk") +
                   b_cost * sq.cost)
  # NL nests & structure
  nlNests = list(root = 1, nSQ = lambdanSQ)
  
  nlStructure = list()
  nlStructure[["root"]] = c("alt3", "nSQ")
  nlStructure[["nSQ"]] = c("alt1", "alt2")
  
  # model components
  nl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V,
    nlNests = nlNests,
    nlStructure = nlStructure
  )
  P[["model"]] = apollo_nl(nl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# estimate model
nl <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# show outputs
apollo_modelOutput(nl,
                   modelOutput_settings=list(printPVal=1))

# write outputs
apollo_saveOutput(nl,
                  saveOutput_settings=list(printPVal=1))

######################
# 5. Multinomial logit with interactions in preference space
######################

apollo_control <- list(
  modelName = "glaciers_CE_mnl_inter",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  outputDirectory = "Model_outputs"
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
                 b_cost = 0,
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

#validate inputs
apollo_inputs <- apollo_validateInputs()

#define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  
  # define interaction terms
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
                   b_climate_1.8 * (sq.climate == "1.8_degrees") +
                   b_flood_freq * (sq.flood == "frequent") +
                   b_drought_freq * (sq.drought == "frequent") +
                   b_water_30 * (sq.freshwater == "30_mil") +
                   b_habitat_high * (sq.habitat == "high_risk") +
                   b_cost * sq.cost)
  
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  )
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# estimate model
mnl_inter <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show outputs
apollo_modelOutput(mnl_inter,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mnl_inter,
                  saveOutput_settings=list(printPVal=1))

######################
# 5a. Multinomial logit with interactions in preference space without insignificant co-variates
######################

apollo_control <- list(
  modelName = "glaciers_CE_mnl_inter_simple",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  outputDirectory = "Model_outputs"
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
                 b_cost = 0,
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
                 # b_male_climate15 = 0,
                 # b_male_climate12 = 0,
                 # b_male_flood10 = 0,
                 b_male_flood20 = 0,
                 # b_male_drought10 = 0,
                 b_male_drought20 = 0,
                 # b_male_water35 = 0,
                 # b_male_water40 = 0,
                 # b_male_habitathalf = 0,
                 # b_urban_climate15 = 0,
                 # b_urban_climate12 = 0,
                 # b_urban_flood10 = 0,
                 # b_urban_flood20 = 0,
                 b_urban_drought10 = 0,
                 # b_urban_drought20 = 0,
                 b_urban_water35 = 0,
                 b_urban_water40 = 0,
                 # b_urban_habitathalf = 0,
                 # b_edu_climate15 = 0,
                 # b_edu_climate12 = 0,
                 b_edu_flood10 = 0,
                 # b_edu_flood20 = 0,
                 b_edu_drought10 = 0,
                 # b_edu_drought20 = 0,
                 b_edu_water35 = 0,
                 b_edu_water40 = 0,
                 # b_edu_habitathalf = 0,
                 # b_income_climate15 = 0,
                 # b_income_climate12 = 0,
                 # b_income_flood10 = 0,
                 # b_income_flood20 = 0,
                 # b_income_drought10 = 0,
                 # b_income_drought20 = 0,
                 # b_income_water35 = 0,
                 # b_income_water40 = 0,
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
  b_flood10 = b_flood_10 + b_basinds_flood10 * (basin == "ds") + b_basinwt_flood10 * (basin == "wt") + 
    b_edu_flood10 * edu + b_exp_flood10 * exp_flood
  b_flood20 = b_flood_20 + b_basinds_flood20 * (basin == "ds") + b_basinwt_flood20 * (basin == "wt") + 
    b_male_flood20 * (gender == 1) + b_exp_flood20 * exp_flood
  b_drought10 = b_drought_10 + b_basinds_drought10 * (basin == "ds") + b_basinwt_drought10 * (basin == "wt") + 
    b_urban_drought10 * (household == 1) + b_edu_drought10 * edu + b_exp_drought10 * exp_drought
  b_drought20 = b_drought_20 + b_basinds_drought20 * (basin == "ds") + b_basinwt_drought20 * (basin == "wt") +
    b_male_drought20 * (gender == 1) + b_exp_drought20 * exp_drought
  b_water35 = b_water_35 + b_basinds_water35 * (basin == "ds") + b_basinwt_water35 * (basin == "wt") + 
    b_urban_water35 * (household == 1) + b_edu_water35 * edu
  b_water40 = b_water_40 + b_basinds_water40 * (basin == "ds") + b_basinwt_water40 * (basin == "wt") +
    b_urban_water40 * (household == 1) + b_edu_water40 * edu
  b_habitathalf = b_habitat_half + b_basinds_habitathalf * (basin == "ds") + b_basinwt_habitathalf * (basin == "wt") + 
    b_income_habitathalf * log(income)
  
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
                   b_climate_1.8 * (sq.climate == "1.8_degrees") +
                   b_flood_freq * (sq.flood == "frequent") +
                   b_drought_freq * (sq.drought == "frequent") +
                   b_water_30 * (sq.freshwater == "30_mil") +
                   b_habitat_high * (sq.habitat == "high_risk") +
                   b_cost * sq.cost)
  
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES,
    V = V
  )
  P[["model"]] = apollo_mnl(mnl_settings, functionality)
  P = apollo_panelProd(P, apollo_inputs, functionality)
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# estimate model
mnl_inter_simple <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

#show outputs
apollo_modelOutput(mnl_inter_simple,
                   modelOutput_settings=list(printPVal=1))

#write outputs
apollo_saveOutput(mnl_inter_simple,
                  saveOutput_settings=list(printPVal=1))

#######################
# 6. Simple mixed logit in WTP space
#######################

## implement a simple mixed logit model to allow for inter-individual preference heterogeneity
apollo_control = list(
  modelName = "glaciers_CE_mxl_wtp",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  mixing = TRUE,
  nCores = 7,
  seed = 2104,
  outputDirectory = "Model_outputs"
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
# 7. Mixed logit with all relevant interactions in preference space
######################

apollo_control <- list(
  modelName = "glaciers_CE_mxl_inter",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  mixing = TRUE,
  nCores = 7,
  seed = 2104,
  outputDirectory = "Model_outputs"
)

# define coefficients to be estimated and starting values for each
# starting values based on simple mixed logit
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
                 b_exp_flood20 = 0,
                 # SQ interactions
                 sq_disappear = 0,
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
  asc_sq = asc_3 + sq_disappear * disappear + sq_way_pollution * way_pollution.control + sq_way_park * way_national.park #+ sq_way_climate
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
  b_water35 = b_water_35 + b_basinds_water35 * (basin == "ds") + b_basinwt_water35 * (basin == "wt") #+
    b_male_water35 * (gender == 1) + b_urban_water35 * (household == 1) + b_edu_water35 * edu + b_income_water35 * income
  b_water40 = b_water_40 + b_basinds_water40 * (basin == "ds") + b_basinwt_water40 * (basin == "wt") #+
    b_male_water40 * (gender == 1) + b_urban_water40 * (household == 1) + b_edu_water40 * edu + b_income_water40 * income
  b_habitathalf = b_habitat_half + b_basinds_habitathalf * (basin == "ds") + b_basinwt_habitathalf * (basin == "wt") #+
    b_male_habitathalf * (gender == 1) + b_urban_habitathalf * (household == 1) + b_edu_habitathalf * edu + b_income_habitathalf * income
  
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
  V[['alt3']] = (asc_sq + 
      b_climate_1.8 * (sq.climate == "1.8_degrees") +
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
# 8. Subregion-specific multinomial logits
######################

# duplicate database
dataset <- database

######
# 8a. Water tower
######
# create database for the subregion
database <- subset(dataset, basin=="wt")

# set required parameters of project
apollo_control <- list(
  modelName = "glaciers_CE_mnl_wt",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  outputDirectory = "Model_outputs"
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

# set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high")

# validate inputs
apollo_inputs <- apollo_validateInputs()

# define model (for details of the required function elements, see apollo documentation)
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
                   b_climate_1.8 * (sq.climate == "1.8_degrees") +
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

# estimate model
mnl_wt <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# show outputs
apollo_modelOutput(mnl_wt,
                   modelOutput_settings=list(printPVal=1))

# write outputs
apollo_saveOutput(mnl_wt,
                  saveOutput_settings=list(printPVal=1))

# calculate WTP mean and standard errors
deltaMethod_settings = list(expression = c(wtp_climate15 = "-b_climate_1.5/b_cost",
                                           wtp_climate12 = "-b_climate_1.2/b_cost",
                                           wtp_flood10 = "-b_flood_10/b_cost",
                                           wtp_flood20 = "-b_flood_20/b_cost",
                                           wtp_drought10 = "-b_drought_10/b_cost",
                                           wtp_drought20 = "-b_drought_20/b_cost",
                                           wtp_water35 = "-b_water_35/b_cost",
                                           wtp_water40 = "-b_water_40/b_cost",
                                           wtp_habhalf = "-b_habitat_half/b_cost"
))
apollo_deltaMethod(mnl_wt, deltaMethod_settings)

######
# 8b. Downstream
######
# create database for the subregion
database <- subset(dataset, basin=="ds")

# set required parameters of project
apollo_control <- list(
  modelName = "glaciers_CE_mnl_ds",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  outputDirectory = "Model_outputs"
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

# set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high")

# validate inputs
apollo_inputs <- apollo_validateInputs()

# define model (for details of the required function elements, see apollo documentation)
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
                   b_climate_1.8 * (sq.climate == "1.8_degrees") +
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

# estimate model
mnl_ds <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# show outputs
apollo_modelOutput(mnl_ds,
                   modelOutput_settings=list(printPVal=1))

# write outputs
apollo_saveOutput(mnl_ds,
                  saveOutput_settings=list(printPVal=1))

# calculate WTP mean and standard errors
deltaMethod_settings = list(expression = c(wtp_climate15 = "-b_climate_1.5/b_cost",
                                           wtp_climate12 = "-b_climate_1.2/b_cost",
                                           wtp_flood10 = "-b_flood_10/b_cost",
                                           wtp_flood20 = "-b_flood_20/b_cost",
                                           wtp_drought10 = "-b_drought_10/b_cost",
                                           wtp_drought20 = "-b_drought_20/b_cost",
                                           wtp_water35 = "-b_water_35/b_cost",
                                           wtp_water40 = "-b_water_40/b_cost",
                                           wtp_habhalf = "-b_habitat_half/b_cost"
))
apollo_deltaMethod(mnl_ds, deltaMethod_settings)

######
# 8c. Other
######
# create database for the subregion
database <- subset(dataset, basin=="oc")

# set required parameters of project
apollo_control <- list(
  modelName = "glaciers_CE_mnl_oc",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  outputDirectory = "Model_outputs"
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

# set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8",
                  "b_flood_freq",
                  "b_drought_freq",
                  "b_water_30",
                  "b_habitat_high")

# validate inputs
apollo_inputs <- apollo_validateInputs()

# define model (for details of the required function elements, see apollo documentation)
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
                   b_climate_1.8 * (sq.climate == "1.8_degrees") +
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

# estimate model
mnl_oc <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# show outputs
apollo_modelOutput(mnl_oc,
                   modelOutput_settings=list(printPVal=1))

# write outputs
apollo_saveOutput(mnl_oc,
                  saveOutput_settings=list(printPVal=1))

# calculate WTP mean and standard errors
deltaMethod_settings = list(expression = c(wtp_climate15 = "-b_climate_1.5/b_cost",
                                           wtp_climate12 = "-b_climate_1.2/b_cost",
                                           wtp_flood10 = "-b_flood_10/b_cost",
                                           wtp_flood20 = "-b_flood_20/b_cost",
                                           wtp_drought10 = "-b_drought_10/b_cost",
                                           wtp_drought20 = "-b_drought_20/b_cost",
                                           wtp_water35 = "-b_water_35/b_cost",
                                           wtp_water40 = "-b_water_40/b_cost",
                                           wtp_habhalf = "-b_habitat_half/b_cost"
))
apollo_deltaMethod(mnl_oc, deltaMethod_settings)

######################
# 9. Latent Class Model
######################

database <- dataset

# set required parameters of project
apollo_control <- list(
  modelName = "glaciers_CE_lcm",
  modelDescr = "Discrete choice experiment into preferences for glacier ES in China",
  indivID = "ID",
  nCores = 7,
  outputDirectory = "Model_outputs"
)

# define coefficients to be estimated and starting values for each
apollo_beta <- c(asc_1 = 0,
                 asc_2 = 0,
                 asc_3 = 0,
                 b_climate_1.8_a = 0,
                 b_climate_1.8_b = 0,
                 b_climate_1.8_c = 0,
                 b_climate_1.5_a = 0,
                 b_climate_1.5_b = 0.05,
                 b_climate_1.5_c = 0.1,
                 b_climate_1.2_a = 0,
                 b_climate_1.2_b = 0.05,
                 b_climate_1.2_c = 0.1,
                 b_flood_freq_a = 0,
                 b_flood_freq_b = 0,
                 b_flood_freq_c = 0,
                 b_flood_10_a = 0,
                 b_flood_10_b = 0.05,
                 b_flood_10_c = 0.1,
                 b_flood_20_a = 0,
                 b_flood_20_b = 0.05,
                 b_flood_20_c = 0.1,
                 b_drought_freq_a = 0,
                 b_drought_freq_b = 0,
                 b_drought_freq_c = 0,
                 b_drought_10_a = 0,
                 b_drought_10_b = 0.05,
                 b_drought_10_c = 0.1,
                 b_drought_20_a = 0,
                 b_drought_20_b = 0.05,
                 b_drought_20_c = 0.1,
                 b_water_30_a = 0,
                 b_water_30_b = 0,
                 b_water_30_c = 0,
                 b_water_35_a = 0,
                 b_water_35_b = 0.05,
                 b_water_35_c = 0.1,
                 b_water_40_a = 0,
                 b_water_40_b = 0.05,
                 b_water_40_c = 0.1,
                 b_habitat_high_a = 0,
                 b_habitat_high_b = 0,
                 b_habitat_high_c = 0,
                 b_habitat_half_a = 0,
                 b_habitat_half_b = 0.05,
                 b_habitat_half_c = 0.1,
                 b_cost_a = 0,
                 b_cost_b = -0.1,
                 b_cost_c = -0.01,
                 # class membership parameters
                 delta_a = 0,
                 gamma_basinds_a = 0,
                 gamma_basinwt_a = 0,
                 gamma_male_a = 0,
                 gamma_urban_a = 0,
                 gamma_edu_a = 0,
                 gamma_income_a = 0,
                 gamma_exp_drought_a = 0,
                 gamma_exp_flood_a = 0,
                 delta_b = 0,
                 gamma_basinds_b = 0,
                 gamma_basinwt_b = 0,
                 gamma_male_b = 0,
                 gamma_urban_b = 0,
                 gamma_edu_b = 0,
                 gamma_income_b = 0,
                 gamma_exp_drought_b = 0,
                 gamma_exp_flood_b = 0,
                 delta_c = 0,
                 gamma_basinds_c = 0,
                 gamma_basinwt_c = 0,
                 gamma_male_c = 0,
                 gamma_urban_c = 0,
                 gamma_edu_c = 0,
                 gamma_income_c = 0,
                 gamma_exp_drought_c = 0,
                 gamma_exp_flood_c = 0
                 )

# set fixed 0 coefficients (here: only one of the alternative-specific constants)
apollo_fixed <- c("asc_1",
                  "b_climate_1.8_a",
                  "b_climate_1.8_b",
                  "b_climate_1.8_c",
                  "b_flood_freq_a",
                  "b_flood_freq_b",
                  "b_flood_freq_c",
                  "b_drought_freq_a",
                  "b_drought_freq_b",
                  "b_drought_freq_c",
                  "b_water_30_a",
                  "b_water_30_b",
                  "b_water_30_c",
                  "b_habitat_high_a",
                  "b_habitat_high_b",
                  "b_habitat_high_c",
                  "delta_a",
                  "gamma_basinds_a",
                  "gamma_basinwt_a",
                  "gamma_male_a",
                  "gamma_urban_a",
                  "gamma_edu_a",
                  "gamma_income_a",
                  "gamma_exp_drought_a",
                  "gamma_exp_flood_a")

# define latent class components
apollo_lcPars = function(apollo_beta, apollo_inputs){
  lcpars = list()
  lcpars[["b_climate_1.8"]] = list(b_climate_1.8_a, b_climate_1.8_b, b_climate_1.8_c)
  lcpars[["b_climate_1.5"]] = list(b_climate_1.5_a, b_climate_1.5_b, b_climate_1.5_c)
  lcpars[["b_climate_1.2"]] = list(b_climate_1.2_a, b_climate_1.2_b, b_climate_1.2_c)
  lcpars[["b_flood_freq"]] = list(b_flood_freq_a, b_flood_freq_b, b_flood_freq_c)
  lcpars[["b_flood_10"]] = list(b_flood_10_a, b_flood_10_b, b_flood_10_c)
  lcpars[["b_flood_20"]] = list(b_flood_20_a, b_flood_20_b, b_flood_20_c)
  lcpars[["b_drought_freq"]] = list(b_drought_freq_a, b_drought_freq_b, b_drought_freq_c)
  lcpars[["b_drought_10"]] = list(b_drought_10_a, b_drought_10_b, b_drought_10_c)
  lcpars[["b_drought_20"]] = list(b_drought_20_a, b_drought_20_b, b_drought_20_c)
  lcpars[["b_water_30"]] = list(b_water_30_a, b_water_30_b, b_water_30_c)
  lcpars[["b_water_40"]] = list(b_water_40_a, b_water_40_b, b_water_40_c)
  lcpars[["b_water_35"]] = list(b_water_35_a, b_water_35_b, b_water_35_c)
  lcpars[["b_habitat_high"]] = list(b_habitat_high_a, b_habitat_high_b, b_habitat_high_c)
  lcpars[["b_habitat_half"]] = list(b_habitat_half_a, b_habitat_half_b, b_habitat_half_c)
  lcpars[["b_cost"]] = list(b_cost_a, b_cost_b, b_cost_c)
  
  # utilities of class allocation model
  V=list()
  V[["class_a"]] = delta_a + gamma_basinds_a * (basin == "ds") + gamma_basinwt_a * (basin == "wt") + 
    gamma_male_a * (gender == 1) + gamma_urban_a * (household == 1) + gamma_edu_a * edu + 
    gamma_income_a * log(income) + gamma_exp_drought_a * exp_drought + gamma_exp_flood_a * exp_flood
  V[["class_b"]] = delta_b + gamma_basinds_b * (basin == "ds") + gamma_basinwt_b * (basin == "wt") + 
    gamma_male_b * (gender == 1) + gamma_urban_b * (household == 1) + gamma_edu_b * edu + 
    gamma_income_b * log(income) + gamma_exp_drought_b * exp_drought + gamma_exp_flood_b * exp_flood
  V[["class_c"]] = delta_c + gamma_basinds_c * (basin == "ds") + gamma_basinwt_c * (basin == "wt") + 
    gamma_male_c * (gender == 1) + gamma_urban_c * (household == 1) + gamma_edu_c * edu + 
    gamma_income_c * log(income) + gamma_exp_drought_c * exp_drought + gamma_exp_flood_c * exp_flood
  
  # settings for class allocation models
  classAlloc_settings = list(
    classes      = c(class_a = 1, class_b = 2, class_c = 3), 
    utilities    = V  
  )
  
  lcpars[["pi_values"]] = apollo_classAlloc(classAlloc_settings)
  
  return(lcpars)
}

# validate inputs
apollo_inputs <- apollo_validateInputs()

# define model (for details of the required function elements, see apollo documentation)
apollo_probabilities <- function(apollo_beta, apollo_inputs, functionality = "estimate"){
  apollo_attach(apollo_beta, apollo_inputs)
  on.exit(apollo_detach(apollo_beta, apollo_inputs))
  P = list()
  mnl_settings = list(
    alternatives = c(alt1 = 1, alt2 = 2, alt3 = 3),
    choiceVar = RES
  )
  # loop over classes
  for (s in 1:3){
    V = list()
    V[['alt1']] = (asc_1 + b_climate_1.8[[s]] * (alt1.climate == "1.8_degrees") + b_climate_1.5[[s]] * (alt1.climate == "1.5_degrees") + 
                     b_climate_1.2[[s]] * (alt1.climate == "1.2_degrees") + b_flood_freq[[s]] * (alt1.flood == "frequent") + 
                     b_flood_10[[s]] * (alt1.flood == "10_%") + b_flood_20[[s]] * (alt1.flood == "20_%") + 
                     b_drought_freq[[s]] * (alt1.drought == "frequent") + b_drought_10[[s]] * (alt1.drought == "10_%") + 
                     b_drought_20[[s]] * (alt1.drought == "20_%") + b_water_30[[s]] * (alt1.freshwater == "30_mil") +
                     b_water_35[[s]] * (alt1.freshwater == "35_mil") + b_water_40[[s]] * (alt1.freshwater == "40_mil") +
                     b_habitat_high[[s]] * (alt1.habitat == "high_risk") + b_habitat_half[[s]] * (alt1.habitat == "half_risk") +
                     b_cost[[s]] * alt1.cost)
    V[['alt2']] = (asc_2 + b_climate_1.8[[s]] * (alt2.climate == "1.8_degrees") + b_climate_1.5[[s]] * (alt2.climate == "1.5_degrees") + 
                     b_climate_1.2[[s]] * (alt2.climate == "1.2_degrees") + b_flood_freq[[s]] * (alt2.flood == "frequent") + 
                     b_flood_10[[s]] * (alt2.flood == "10_%") + b_flood_20[[s]] * (alt2.flood == "20_%") + 
                     b_drought_freq[[s]] * (alt2.drought == "frequent") + b_drought_10[[s]] * (alt2.drought == "10_%") + 
                     b_drought_20[[s]] * (alt2.drought == "20_%") + b_water_30[[s]] * (alt2.freshwater == "30_mil") +
                     b_water_35[[s]] * (alt2.freshwater == "35_mil") + b_water_40[[s]] * (alt2.freshwater == "40_mil") +
                     b_habitat_high[[s]] * (alt2.habitat == "high_risk") + b_habitat_half[[s]] * (alt2.habitat == "half_risk") +
                     b_cost[[s]] * alt2.cost)
    V[['alt3']] = (asc_3 + b_climate_1.8[[s]] * (sq.climate == "1.8_degrees") + b_climate_1.5[[s]] * (sq.climate == "1.5_degrees") + 
                     b_climate_1.2[[s]] * (sq.climate == "1.2_degrees") + b_flood_freq[[s]] * (sq.flood == "frequent") + 
                     b_flood_10[[s]] * (sq.flood == "10_%") + b_flood_20[[s]] * (sq.flood == "20_%") + 
                     b_drought_freq[[s]] * (sq.drought == "frequent") + b_drought_10[[s]] * (sq.drought == "10_%") + 
                     b_drought_20[[s]] * (sq.drought == "20_%") + b_water_30[[s]] * (sq.freshwater == "30_mil") +
                     b_water_35[[s]] * (sq.freshwater == "35_mil") + b_water_40[[s]] * (sq.freshwater == "40_mil") +
                     b_habitat_high[[s]] * (sq.habitat == "high_risk") + b_habitat_half[[s]] * (sq.habitat == "half_risk") +
                     b_cost[[s]] * sq.cost)
    
    mnl_settings$utilities = V
    
    # compute within-class choice probabilities using MNL model
    P[[paste0("Class_",s)]] = apollo_mnl(mnl_settings, functionality)
    
    # take product across observation for same individual
    P[[paste0("Class_",s)]] = apollo_panelProd(P[[paste0("Class_",s)]], apollo_inputs ,functionality)
  }
  
  # compute latent class model probabilities
  lc_settings  = list(inClassProb = P, classProb = pi_values)
  P[["model"]] = apollo_lc(lc_settings, apollo_inputs, functionality)
  
  # prepare and return outputs of function
  P = apollo_prepareProb(P, apollo_inputs, functionality)
  return(P)
}

# optional
apollo_beta <- apollo_searchStart(apollo_beta, apollo_fixed,apollo_probabilities, apollo_inputs)

# estimate model
lcm <- apollo_estimate(apollo_beta, apollo_fixed, apollo_probabilities, apollo_inputs)

# show outputs
apollo_modelOutput(lcm,
                   modelOutput_settings=list(printPVal=1))

# write outputs
apollo_saveOutput(lcm,
                  saveOutput_settings=list(printPVal=1))

######################
# 10. Save outputs as DOCX tables
######################

## Main results
# optionally: load required models
mnl_inter_simple <- apollo_loadModel("glaciers_CE_mnl_inter_simple")
mnl <- apollo_loadModel("glaciers_CE_mnl")

# save model as texreg-readable object
mnl_inter_simple2 <- quicktexregapollo(mnl_inter_simple)
mnl2 <- quicktexregapollo(mnl)

# write output to word table
wordreg(list(mnl2, mnl_inter_simple2), "Tables/glaciers_CE_main.docx",
        digits = 3, single.row = T)

## Supplementary material
# Multinomial logit with all interactions

# mnl_inter <- apollo_loadModel("glaciers_CE_mnl_inter")

mnl_inter2 <- quicktexregapollo(mnl_inter)

wordreg(mnl_inter2, "Tables/glaciers_CE_mnl_inter.docx",
        digits = 3, single.row = T)

# Multinomial logit for subregions

# mnl_wt <- apollo_loadModel("glaciers_CE_mnl_wt")
# mnl_ds <- apollo_loadModel("glaciers_CE_mnl_ds")
# mnl_oc <- apollo_loadModel("glaciers_CE_mnl_oc")

mnl.wt <- quicktexregapollo(mnl_wt)
mnl.ds <- quicktexregapollo(mnl_ds)
mnl.oc <- quicktexregapollo(mnl_oc)

wordreg(list(mnl.wt, mnl.ds, mnl.oc), "Tables/glaciers_CE_mnl_basins.docx",
        digits = 3, single.row = T)

# Nested logit

# nl <- apollo_loadModel("glaciers_CE_nl")

nl2 <- quicktexregapollo(nl)

wordreg(nl2, "Tables/glaciers_CE_nl.docx",
        digits = 3, single.row = T)

# Mixed logits (simple and with interactions)

# mxl <- apollo_loadModel("glaciers_CE_mxl")
# mxl_inter <- apollo_loadModel("glaciers_CE_mxl_inter")

mxl2 <- quicktexregapollo(mxl)
mxl_inter2 <- quicktexregapollo(mxl_inter)

wordreg(list(mxl2, mxl_inter2), "Tables/glaciers_CE_mxl.docx",
        digits = 3, single.row = T)
