# ----------- #
# OBJECTIVE: Estimation of the ODE-based model using viral load data (without immuno markers), with or without group effects
# Author: Marie Alexandre
# Date: 2021/11/23

# R version: 3.6.2
# Monolix version: 2019R1

# --- Description: 
# In this file, we estimate via Monolix the ODE-based model describing the viral load dynamics (gRNA and sgRNA) within both trachea and nasopharynx comprtments
# Data simulated using the R file "Rcode/Data_Simulation.R" are used in the file (Immuno markers are not used in this part)
# Immuno markers being not considered in this file, group effect are assumed on both parameters beta and delta
# In this file, 4 models are estimated:
#   a) Model without group covariates
#   b) Model with group covariate on beta
#   c) Model with group covariate on delta
#   d) Model with group covariates on both beta and delta

# Tips: Ctrl + shift + O to extend document outline
# ----------- #


rm(list=ls())
'%notin%' <- Negate('%in%')


# Libraries --- ####
require(Rsmlx) ; require(lixoftConnectors)
require(plyr)
# ----- #

# Initialization of the lixoft connectors (for simulx) --- ####
Monolix_path <- "C:/ProgramData/Lixoft/MonolixSuite2019R1"  # To modify if necessary
initializeLixoftConnectors(software = "monolix", path = Monolix_path,force=TRUE)
# ----- #


Data_Folder <- "Simulated_data"
data_file <- "Simulated_data_ViralLoad_ImmunoMarkers" # Only viral load are used in this file
Model_Folder <- "Mlxtran_Models"
Model_file <- paste(Model_Folder,"MlxtranCode_MechanisticModel_SARSCoV2_VL_NoExternalCov.txt",sep="/")

Monolix_Folder <- "Monolix_Estimation"

# Model Specifications
Initial_parameters <- data.frame(name=c("beta_N_pow","fact_beta_T","fact_delta_T","delta_N","c","cI","k","g","mu","P_N","fact_P_T","alpha_VLSG","thresh_Weight","omega_beta_N_pow","omega_delta_N"),
                                 values=c(-7,0,0,1,3,20,3,0,0.001,40000,-1,1.33,4.5,1,1),
                                 method=c("MLE","FIXED","FIXED","MLE","FIXED","FIXED","FIXED","FIXED","FIXED","MLE","MLE","FIXED","FIXED","MLE","MLE"),
                                 stringsAsFactors = FALSE)

# Definition of parameter variability
parameter_variability <- c(TRUE,FALSE,FALSE,TRUE,rep(FALSE,9))
names(parameter_variability) <- c("beta_N_pow","fact_beta_T","fact_delta_T","delta_N","c","cI","k","g","mu","P_N","fact_P_T","alpha_VLSG","thresh_Weight")

# Definition individual parameters distribution 
parameter_distribution <- c("normal","normal","normal","logNormal","logNormal","logNormal","logNormal","normal","lognormal","logNormal","normal","lognormal","lognormal")
names(parameter_distribution) <-  c("beta_N_pow","fact_beta_T","fact_delta_T","delta_N","c","cI","k","g","mu","P_N","fact_P_T","alpha_VLSG","thresh_Weight")

# Definition of the header of dataset in Monolix project
Header_Monolix <- c("id","time","observation","obsid","catcov","ignore","ignore","regressor","cens","ignore","ignore","ignore",rep("ignore",6))
names(Header_Monolix) <- c("AnimalID","Time","log10VL","obs","Group","type","organ","Weight","censored","ECL","ECL.Slope","ECL.intercept","Marker2","Marker2.Slope","Marker2.Intercept","Marker3","Marker3.Slope","Marker3.Intercept")

# Definition of the type of data
obs_type <- rep("continuous",4) 
names(obs_type) <- c("y_1","y_2","y_3","y_4")
observationNames <- c("y_1","y_2","y_3","y_4")



# ------------------------- #
# PART I: Estimation of the model without any covariates (Group or markers) ####
Project_Name_NoCov <- "Model_Estimation_NoCov"


# Creation of the Monolix project
newProject(data=list(dataFile=paste(Data_Folder,paste(data_file,".csv",sep=""),sep="/"),
                     header=names(Header_Monolix),
                     headerTypes=as.vector(Header_Monolix),
                     observationNames=observationNames,
                     observationTypes = obs_type,
                     mapping=list('1'="y_1",'2'="y_2",'3'="y_3",'4'="y_4")),
           modelFile=Model_file)
saveProject(projectFile = paste(Monolix_Folder,paste(Project_Name_NoCov,".mlxtran",sep=""),sep="/"))


# Modification of the project settings
Scenario <- getScenario()
Scenario$tasks <- c(populationParameterEstimation=T,conditionalDistributionSampling=T,
                    conditionalModeEstimation=T,standardErrorEstimation=T,
                    logLikelihoodEstimation=T,plots=FALSE)
setScenario(Scenario)

# Modification of the model specification
setErrorModel(y1="constant",y2="constant",y3="constant",y4="constant")
setIndividualParameterVariability(parameter_variability)
setIndividualParameterDistribution(parameter_distribution)

# Modification of the initial condition
Init_params <- getPopulationParameterInformation()
for(p in 1:nrow(Init_params)){
  param <- Init_params$name[p]
  selected_param <- strsplit(param,split="_pop",fixed=TRUE)[[1]]
  if(selected_param %in% Initial_parameters$name){
    ind <- which(Initial_parameters$name == selected_param)
    Init_params$initialValue[p] <- Initial_parameters$values[ind]
    Init_params$method[p] <- Initial_parameters$method[ind]
  }
} 
setPopulationParameterInformation(Init_params)
setPopulationParameterEstimationSettings(exploratoryautostop=FALSE,nbexploratoryiterations=800)
setPopulationParameterEstimationSettings(smoothingautostop=FALSE,nbsmoothingiterations=200)

setStandardErrorEstimationSettings(minIterations=100,maxIterations=5000)
saveProject()

# Model Estimation
runScenario()
saveProject()
# ------------------------- #


# ------------------------- #
# PART II: Estimation of the model with group effect on beta only ####
Project_Name_GroupBeta <- "Model_Estimation_CovGroupBeta"

# We load the model without covariate and we add the covariates "Group" on beta and delta
loadProject(projectFile = paste(Monolix_Folder,paste(Project_Name_NoCov,".mlxtran",sep=""),sep="/"))
saveProject(projectFile = paste(Monolix_Folder,paste(Project_Name_GroupBeta,".mlxtran",sep=""),sep="/"))

# Addition of the group covariates
addCategoricalTransformedCovariate(groupVsnaive=list(reference="G_Naive",from="Group",transformed=list(G_Naive="Naive",G_conv="Convalescent",G_convCD40="Conv-CD40")))
setCovariateModel(beta_N_pow=c(groupVsnaive=TRUE))
saveProject()

# Model Estimation
runScenario()
saveProject()
# ------------------------- #

# ------------------------- #
# PART III: Estimation of the model with group effect on delta only ####
Project_Name_GroupDelta <- "Model_Estimation_CovGroupDelta"

# We load the model without covariate and we add the covariates "Group" on beta and delta
loadProject(projectFile = paste(Monolix_Folder,paste(Project_Name_NoCov,".mlxtran",sep=""),sep="/"))
saveProject(projectFile = paste(Monolix_Folder,paste(Project_Name_GroupDelta,".mlxtran",sep=""),sep="/"))

# Addition of the group covariates
addCategoricalTransformedCovariate(groupVsnaive=list(reference="G_Naive",from="Group",transformed=list(G_Naive="Naive",G_conv="Convalescent",G_convCD40="Conv-CD40")))
setCovariateModel(delta_N=c(groupVsnaive=TRUE))
saveProject()

# Model Estimation
runScenario()
saveProject()
# ------------------------- #

# ------------------------- #
# PART IV: Estimation of the model with group effects on both beta and delta  ####
Project_Name_GroupBetaDelta <- "Model_Estimation_CovGroupBetaDelta"

# We load the model without covariate and we add the covariates "Group" on beta and delta
loadProject(projectFile = paste(Monolix_Folder,paste(Project_Name_NoCov,".mlxtran",sep=""),sep="/"))
saveProject(projectFile = paste(Monolix_Folder,paste(Project_Name_GroupBetaDelta,".mlxtran",sep=""),sep="/"))

# Addition of the group covariates
addCategoricalTransformedCovariate(groupVsnaive=list(reference="G_Naive",from="Group",transformed=list(G_Naive="Naive",G_conv="Convalescent",G_convCD40="Conv-CD40")))
setCovariateModel(beta_N_pow=c(groupVsnaive=TRUE),delta_N=c(groupVsnaive=TRUE))
saveProject()

# Model Estimation
runScenario()
saveProject()
# ------------------------- #