# --------------- # 
# DESCRIPTION: Estimation of the mechanistic Joint model describing both viral and antibody dynamics following infection
#
# Author: Marie Alexandre 
# R version: 4.2.1
# Monolix Version: 2023R1


# Open the document outline (Ctrl + shift + O) to see the structure of the document
# --------------- # 





rm(list=ls())
'%notin%' <- Negate('%in%')




# --- LIBRAIRIES ----
library(lixoftConnectors,lib.loc = "C:/Users/marie/AppData/Local/R/win-library/4.4/Monolix2023")
# ------------------ #


# --- FUNCTIONS ----
Modification_Initial_condition <- function(currentValues,newValues){
  if("MAP" %in% newValues$method){
    for(p in 1:nrow(currentValues)){
      # p <- 1
      param <- currentValues$name[p]
      selected_param <- strsplit(param,split="_pop",fixed=TRUE)[[1]]
      if(selected_param %in% newValues$name){
        ind <- which(newValues$name == selected_param)
        currentValues$initialValue[p] <- newValues$values[ind]
        currentValues$method[p] <- newValues$method[ind]
        currentValues$priorValue[p] <- newValues$priorValue[ind]
        currentValues$priorSD[p] <- newValues$priorSD[ind]
      }
    }
  }else{
    for(p in 1:nrow(currentValues)){
      # p <- 1
      param <- currentValues$name[p]
      selected_param <- strsplit(param,split="_pop",fixed=TRUE)[[1]]
      if(selected_param %in% newValues$name){
        ind <- which(newValues$name == selected_param)
        currentValues$initialValue[p] <- newValues$values[ind]
        currentValues$method[p] <- newValues$method[ind]
      }
    }
  }
  return(currentValues)
}
# ------------------ #


# --- General variables ----
Project_Folder <- "EarlyViralandAntibodyDynPostInfection"

Data_Folder <- "Simulated_data"
Data_file <- "Simulated_dataset_VL_BAb_NAb_Delta_PostExpo.txt"

Monolix_Folder <- "Model_Estimation_Monolix"
Mlxtran_file <- "MlxtranCode_NHP_JointModel_VL_Ab_Monolix.txt"

MonolixProject <- "Final_JointModel_NHP_DeltaInfection_MONOLIX.mlxtran"
# ---------------- #




Observed_data <- read.table(file = paste(Project_Folder,Data_Folder,Data_file,sep="/"),header=TRUE,sep="\t",dec=".")


# --- Creation of Monolix project ----
initializeLixoftConnectors(software = "monolix",path = "C:/ProgramData/Lixoft/MonolixSuite2023R1",force=TRUE)

# Definition of data headers
Data_header <- c(SubjectID = "id",Group = "catcov",StatutImm = "ignore",DayPostExpo = "time",
                 Value = "ignore",Censored = "cens",CensoredValue = "ignore", MeasureType = "ignore",
                 SampleType = "ignore", Observation = "ignore",CensObservation = "observation", obsid = "obsid",
                 Inoc = "regressor", Inoc_Voc= "ignore",Weight = "regressor", BAb0 = "ignore",NAb0 = "ignore")

# Creation of the Monolix project
newProject(data=list(dataFile = paste(Project_Folder,Data_Folder,Data_file,sep="/"),
                     headerTypes = as.character(Data_header),
                     observationTypes = list("1" = "continuous", "2" = "continuous", 
                                             "3" = "continuous", "4" = "continuous",
                                             "5" = "continuous", "6" = "continuous"),
                     mapping = list(list(data = "1", prediction = "gRNA_T", model = "y1"),
                                    list(data = "2", prediction = "sgRNA_T", model = "y2"),
                                    list(data = "3", prediction = "gRNA_N", model = "y3"),
                                    list(data = "4", prediction = "sgRNA_N", model = "y4"),
                                    list(data = "5", prediction = "Obs_BAb", model = "y5"),
                                    list(data = "6", prediction = "Obs_NAb", model = "y6"))),
           modelFile=paste(Project_Folder,Monolix_Folder,Mlxtran_file,sep="/"))

saveProject(projectFile = paste(Project_Folder,Monolix_Folder,MonolixProject,sep="/"))
# ---------------- #




# --- Modification of estimation features ----
## 1. Model specification  ----

###  > Error model ----
setErrorModel(y1="constant",y2="constant",y3="constant",y4="constant",y5="constant",y6="constant")
setObservationDistribution(y1="normal",y2="normal",y3="normal",y4="normal",y5="normal",y6="normal")


### > Model parameters distribution ----
# Initially, all parameters are used as logNormal. We only modify here parameters that are not logNormal
Parameter_Distribution <- c(fact_beta_T="normal",fact_delta_T="normal",fact_P_T="normal")
setIndividualParameterDistribution(Parameter_Distribution)

###  > Model parameters variability ----
# Initially all parameters are adjusted for inter-individual variability.
Parameter_variability <- getIndividualParameterModel()$variability$id
# We first remove variability from all parameters
Parameter_variability <- setNames(as.logical(1 - Parameter_variability),names(Parameter_variability))
# Addition of variability on delta, S0, theta, eta
Parameter_variability[names(Parameter_variability) %in% c("delta_N","S0","theta_BAb","eta")] <- TRUE
setIndividualParameterVariability(Parameter_variability)

### > Definition of new covariates ----
GroupEffect = list(reference = "Naive", from = "Group", 
                   transformed = list(Naive = c("Naive"), CD40RBDv = c("CD40.RBDv - Naive"),
                                      Conv = c("Conv"), CD40RBDvConv = c("CD40.RBDv - Conv"),
                                      mRNAConv = c("mRNA - Conv"), PanCoVConv = c("CD40.PanCoV - Conv")))

CD40RBDvNaive = list(reference = "0", from = "Group", 
                     transformed = list(`0` = c("Naive","Conv","CD40.RBDv - Conv","CD40.PanCoV - Conv","mRNA - Conv"),
                                        `1` = c("CD40.RBDv - Naive")))

Conv = list(reference = "0", from = "Group", 
            transformed = list(`0` = c("Naive","CD40.RBDv - Naive","CD40.RBDv - Conv","CD40.PanCoV - Conv","mRNA - Conv"),
                               `1` = c("Conv")))

CD40RBDvConv = list(reference = "0", from = "Group", 
                    transformed = list(`0` = c("Naive","CD40.RBDv - Naive","Conv","CD40.PanCoV - Conv","mRNA - Conv"),
                                       `1` = c("CD40.RBDv - Conv")))

CD40PanCovConv = list(reference = "0", from = "Group", 
                      transformed = list(`0` = c("Naive","CD40.RBDv - Naive","Conv","CD40.RBDv - Conv","mRNA - Conv"),
                                         `1` = c("CD40.PanCoV - Conv")))

mRNAConv = list(reference = "0", from = "Group", 
                transformed = list(`0` = c("Naive","CD40.RBDv - Naive","Conv","CD40.RBDv - Conv","CD40.PanCoV - Conv"),
                                   `1` = c("mRNA - Conv")))

Convalescence = list(reference = "0", from = "Group",
                     transformed = list(`0` = c("Naive","CD40.RBDv - Naive"),
                                        `1` = c("Conv","CD40.RBDv - Conv","CD40.PanCoV - Conv","mRNA - Conv")))

ConvXVacc = list(reference = "0", from ="Group",
                 transformed = list(`0` = c("Naive","CD40.RBDv - Naive","Conv"),
                                    CD40 = c("CD40.RBDv - Conv"),
                                    PanCov = c("CD40.PanCoV - Conv"),
                                    mRNA = c("mRNA - Conv")))

addCategoricalTransformedCovariate(GroupEffect = GroupEffect, CD40RBDvNaive = CD40RBDvNaive,
                                   Conv = Conv, CD40RBDvConv = CD40RBDvConv,
                                   CD40PanCovConv = CD40PanCovConv, mRNAConv = mRNAConv,
                                   Convalescence = Convalescence, ConvXVacc = ConvXVacc)

setCovariateModel(delta_N = c(Conv = TRUE), S0 = c(GroupEffect = TRUE), rho = c(Convalescence = TRUE), 
                  theta_BAb = c(Conv = TRUE, CD40RBDvNaive = TRUE),eta = c(ConvXVacc = TRUE))




## 2. Modification of initial conditions  ----
Current_Init_params <- getPopulationParameterInformation()

# Modification of values
Chosen_initial_values <- rbind(data.frame(name="betas_N",values=10,method="MLE",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="fact_beta_T",values=0,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="delta_N",values=1.0,method="MLE",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="fact_delta_T",values=0,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="P_N",values=10000,method="MLE",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="fact_P_T",values=0,method="MLE",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="alpha_VLSG",values=1.4,method="MLE",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="c",values=3,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="cI",values=20,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="k",values=3,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="mu",values=1E-3,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               
                               data.frame(name="S0",values=1,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F), # Value for naive NHPs
                               data.frame(name="beta_S0_GroupEffect_Conv",values=0,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F), # No effect of Conv on S0
                               data.frame(name="gamma",values=1E-5,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="rho",values=1,method="MLE",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="Smax",values=1000,method="MLE",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="tauS",values=2.3,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="theta_BAb",values=5,method="MLE",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="alpha_NAb",values=1E-4,method="MLE",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="tauAb",values=12,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               data.frame(name="eta",values=1.0,method="MLE",priorValue=NA,priorSD=NA,stringsAsFactors=F),
                               
                               data.frame(name="thresh_Weight",values=4.5,method="FIXED",priorValue=NA,priorSD=NA,stringsAsFactors=F))

Current_Init_params <- Modification_Initial_condition(currentValues=Current_Init_params,newValues=Chosen_initial_values)
setPopulationParameterInformation(Current_Init_params)


## 3. Modification  of estimation settings ----
setPopulationParameterEstimationSettings(exploratoryautostop=TRUE,nbexploratoryiterations=1000,variability="none")
setStandardErrorEstimationSettings(minIterations=100,maxIterations=5000)
setProjectSettings(seed = 123456)
setGeneralSettings(nbChains = 5)


## 4. Definition of estimation scenario ----
Scenario <- getScenario()
Scenario$tasks <- c(populationParameterEstimation=T,conditionalDistributionSampling=T,
                    conditionalModeEstimation=T,standardErrorEstimation=T,
                    logLikelihoodEstimation=T,plots=F)
Scenario$linearization <- FALSE
setScenario(Scenario)

saveProject()
# ---------------- #


# --- Run model estimation ----
runScenario()
# ---------------- #
