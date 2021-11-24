# ----------- #
# OBJECTIVE: Definition of the Monolix model specification for the first step of the covariate selection
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
#   e. the dataset header used in Monolix according to already fixed covariates
#   f. the model (txt files) to use for testing adjustment of covariates on each parameter 

# The parameter beta_immuno correspond to the regression coefficient involved in the adjustment of the parameter of interest for the tested marker
# ----------- #

# a. Definition of the initial parameters
Initial_parameters <- data.frame(name=c("beta_N_pow","fact_beta_T","fact_delta_T","delta_N","c","cI","k","g","mu","P_N","fact_P_T","alpha_VLSG","thresh_Weight","beta_immuno","omega_beta_N_pow","omega_delta_N"),
                                 values=c(-7,0,0,1.24154,3,20,3,0,0.001,40000,0,1,4.5,0,1,1),
                                 method=c("MLE","FIXED","FIXED","MLE","FIXED","FIXED","FIXED","FIXED","FIXED","MLE","MLE","MLE","FIXED","MLE","MLE","MLE"),
                                 stringsAsFactors = FALSE)

# b. Definition of parameter variability
parameter_variability <- c(TRUE,FALSE,FALSE,TRUE,rep(FALSE,10))
names(parameter_variability) <- c("beta_N_pow","fact_beta_T","fact_delta_T","delta_N","c","cI","k","g","mu","P_N","fact_P_T","alpha_VLSG","thresh_Weight","beta_immuno")


# c. Definition individual parameters distribution 
parameter_distribution <- c("normal","normal","normal","logNormal","logNormal","logNormal","logNormal","normal","lognormal","logNormal","normal","lognormal","lognormal","normal")
names(parameter_distribution) <-  c("beta_N_pow","fact_beta_T","fact_delta_T","delta_N","c","cI","k","g","mu","P_N","fact_P_T","alpha_VLSG","thresh_Weight","beta_immuno")


# d. list of covariates (among immuno markers) already added in the model
# At the first step of the covariate selection, we mostly consider no adjustment of model parameters for markers
Fixed_timeCov <- list(beta=NULL,delta=NULL,P=NULL,c=NULL)

# e. Definition of the Monolix Header
if(is.null(unlist(Fixed_timeCov))){
  Header_Monolix <- c("time","catcov","id","regressor","obsid","observation","cens")
  names(Header_Monolix) <- c("Time","Group","AnimalID","Weight","obs","log10VL","censored")
}else{
  Header_Monolix <- c("time","catcov","id","regressor","obsid","observation","cens",rep("regressor",2*length(unlist(Fixed_timeCov))))
  names(Header_Monolix) <- c("Time","Group","AnimalID","Weight","obs","log10VL","censored",colnames(New_dataSet)[8:(7+2*length(unlist(Fixed_timeCov)))])
}

# f. list of models
Model_Files <- list(beta="MlxtranCode_MechanisticModel_SARSCoV2_VL_TimeVaryingCovBeta.txt",
                    delta="MlxtranCode_MechanisticModel_SARSCoV2_VL_TimeVaryingCovDelta.txt")
