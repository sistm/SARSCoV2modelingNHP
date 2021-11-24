# ----------- #
# OBJECTIVE: Definition of the Monolix model specification for the second step of the covariate selection
# Author: Marie Alexandre
# Date: 2021/11/23

# R version: 3.6.2
# Monolix version: 2019R1

# --- Description: 
# In this file, we define:
#   a. the initial values of the model parameters for the model estimation
#   b. the inter-individual variability
#   c. the distribution of the parameters
#   d. the initial convariate structure (among immuno markers)
#   e. the model (txt files) to use for testing adjustment of covariates on each parameter 

# The parameters beta_immuno1 and beta_immuno2 correspond to the regression coefficients involved in the adjustment selected at the first step and of the parameter of interest for the tested marker, respectively
# ----------- #


# a. Definition of the initial parameters
Initial_parameters <- data.frame(name=c("beta_N_pow","fact_beta_T","fact_delta_T","delta_N","c","cI","k","g","mu","P_N","fact_P_T","alpha_VLSG","thresh_Weight","beta_immuno1","beta_immuno2","omega_beta_N_pow","omega_delta_N"),
                                 values=c(-7,0,0,1.24154,3,20,3,0,0.001,40000,0,1,4.5,0,0,1,1),
                                 method=c("MLE","FIXED","FIXED","MLE","FIXED","FIXED","FIXED","FIXED","FIXED","MLE","MLE","MLE","FIXED","MLE","MLE","MLE","MLE"),
                                 stringsAsFactors = FALSE)

# b. Definition of parameter variability
parameter_variability <- c(TRUE,FALSE,FALSE,TRUE,rep(FALSE,11))
names(parameter_variability) <- c("beta_N_pow","fact_beta_T","fact_delta_T","delta_N","c","cI","k","g","mu","P_N","fact_P_T","alpha_VLSG","thresh_Weight","beta_immuno1","beta_immuno2")

# c. Definition individual parameters distribution 
parameter_distribution <- c("normal","normal","normal","logNormal","logNormal","logNormal","logNormal","normal","lognormal","logNormal","normal","lognormal","lognormal","normal","normal")
names(parameter_distribution) <-  c("beta_N_pow","fact_beta_T","fact_delta_T","delta_N","c","cI","k","g","mu","P_N","fact_P_T","alpha_VLSG","thresh_Weight","beta_immuno1","beta_immuno2")

# d. list of covariates (among immuno markers) already added in the model
# The parameter is now adjusted for ECL
Fixed_timeCov <- list(beta=c("ECL"),delta=NULL,P=NULL,c=NULL)

# e. list of models
Model_Files <- list(beta="MlxtranCode_MechanisticModel_SARSCoV2_VL_2TimeVaryingCovBeta.txt",
                    delta="MlxtranCode_MechanisticModel_SARSCoV2_VL_1TimeVaryingCovBeta_1TimeVaryingCovDelta.txt")