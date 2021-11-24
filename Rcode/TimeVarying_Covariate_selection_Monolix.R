# ----------- #
# OBJECTIVE: Stepwise Selection of Time-varying covariates (Immuno markers)
# Author: Marie Alexandre
# Date: 2021/11/23

# R version: 3.6.2
# Monolix version: 2019R1

# --- Description: 
# In this file, we implemented the covariate selection approach developed to identify mechanistic correlate of protection 
# Different type of data transformation are proposed (in addition to the use of natural scale data)
# Proposed transformations: Min-Max,log10,BoxCox (EnvStats package),YeoJohnson (bestNormalize package),Arcsinh (bestNormalize package),Log (bestNormalize package)

# Tips: Ctrl + shift + O to extend document outline
# ----------- #



rm(list=ls())
'%notin%' <- Negate('%in%')


# Libraries --- ####
require(MASS) ; require(reshape2) ; require(dplyr)
require(Rsmlx) ; require(lixoftConnectors)
require(parallel) ; require(snow) ; require(doSNOW)
# ----- #

# Initialization of the lixoft connectors (for simulx) --- ####
Monolix_path <- "C:/ProgramData/Lixoft/MonolixSuite2019R1"  # To modify if necessary
initializeLixoftConnectors(software = "monolix", path = Monolix_path,force=TRUE)
# ----- #

# Functions --- ####
# Data transformation
Data_transformation_function <- function(Data,markers,Transf_name){
  ind_col <- which(colnames(Data) %in% markers)
  Reshaped_data <- Data
  estimated_params <- NULL
  
  for(c in ind_col){
    if(Transf_name == "MinMax_transf"){
      Reshaped_data[,c] <- MinMax_transformation(Data[,c]) 
    }else if(Transf_name == "Log10"){
      Reshaped_data[,c] <- Log10_transformation(Data[,c])
    }else if(Transf_name == "BoxCox"){
      tmp_data <- unique(Reshaped_data[,c("AnimalID","Time",names(Reshaped_data)[c])])
      tmp_data <- tmp_data[which(!is.na(tmp_data[,3]) & tmp_data$Time>=0),]
      boxcox_res <- BoxCox_transformation(values=tmp_data[,3],optimize=TRUE)
      estimated_params <- c(estimated_params,boxcox_res$lambda)
      tmp_data$transf_data <- boxcox_res$data
    }else if(Transf_name == "YeoJohnson-Norm"){
      tmp_data <- unique(Reshaped_data[,c("AnimalID","Time",names(Reshaped_data)[c])])
      tmp_data <- tmp_data[which(!is.na(tmp_data[,3]) & tmp_data$Time>=0),]
      yeojohnson_res <- YeoJohnson_transformation(values=tmp_data[,3],standard=TRUE)
      estimated_params <- c(estimated_params,yeojohnson_res$lambda)
      tmp_data$transf_data <- yeojohnson_res$data
    }else if(Transf_name == "YeoJohnson"){
      tmp_data <- unique(Reshaped_data[,c("AnimalID","Time",names(Reshaped_data)[c])])
      tmp_data <- tmp_data[which(!is.na(tmp_data[,3]) & tmp_data$Time>=0),]
      yeojohnson_res <- YeoJohnson_transformation(values=tmp_data[,3],standard=FALSE)
      estimated_params <- c(estimated_params,yeojohnson_res$lambda)
      tmp_data$transf_data <- yeojohnson_res$data
    }else if(Transf_name == "Arcsinh"){
      tmp_data <- unique(Reshaped_data[,c("AnimalID","Time",names(Reshaped_data)[c])])
      tmp_data <- tmp_data[which(!is.na(tmp_data[,3]) & tmp_data$Time>=0),]
      arcsinh_res <- Arcsinh_transformation(values=tmp_data[,3],standard=FALSE)
      tmp_data$transf_data <- arcsinh_res$data
    }else if(Transf_name == "Arcsinh-Norm"){
      tmp_data <- unique(Reshaped_data[,c("AnimalID","Time",names(Reshaped_data)[c])])
      tmp_data <- tmp_data[which(!is.na(tmp_data[,3]) & tmp_data$Time>=0),]
      arcsinh_res <- Arcsinh_transformation(values=tmp_data[,3],standard=TRUE)
      tmp_data$transf_data <- arcsinh_res$data
    }else if(Transf_name == "Sqrt"){
      tmp_data <- unique(Reshaped_data[,c("AnimalID","Time",names(Reshaped_data)[c])])
      tmp_data <- tmp_data[which(!is.na(tmp_data[,3]) & tmp_data$Time>=0),]
      sqrt_res <- sqrt(tmp_data[,3])
      tmp_data$transf_data <- sqrt_res
    }else if(Transf_name == "Fourthroot"){
      tmp_data <- unique(Reshaped_data[,c("AnimalID","Time",names(Reshaped_data)[c])])
      tmp_data <- tmp_data[which(!is.na(tmp_data[,3]) & tmp_data$Time>=0),]
      fourthroot_res <- (tmp_data[,3])^(0.25)
      tmp_data$transf_data <- fourthroot_res
    }else if(Transf_name == "Log"){
      tmp_data <- unique(Reshaped_data[,c("AnimalID","Time",names(Reshaped_data)[c])])
      tmp_data <- tmp_data[which(!is.na(tmp_data[,3]) & tmp_data$Time>=0),]
      log_res <- Log_transformation(values=tmp_data[,3],standard=FALSE)
      estimated_params <- c(estimated_params,log_res$a)
      tmp_data$transf_data <- log_res$data
    }
    
    if(Transf_name %in% c("BoxCox","YeoJohnson-Norm","YeoJohnson","Log")){
      Reshaped_data[,c] <- sapply(seq(1,nrow(Reshaped_data)),function(i){
        # i <- 1
        tmp_res <- tmp_data$transf_data[which(tmp_data$AnimalID == Reshaped_data$AnimalID[i] & tmp_data$Time == Reshaped_data$Time[i])]
        if(length(tmp_res) == 0){
          tmp_res <- NA
        }
        return(tmp_res)
      })
    }
  }
  
  if(Transf_name %in% c("BoxCox","YeoJohnson-Norm","YeoJohnson","Log")){
    names(estimated_params) <- markers
    return(list(data=Reshaped_data,lambda=estimated_params))
  }else{
    return(Reshaped_data)
  }
}
MinMax_transformation <- function(values){
  res <- (values-min(values,na.rm=TRUE))/(max(values,na.rm=TRUE)-min(values,na.rm=TRUE))
  return(res)
}
Log10_transformation <- function(values){
  res <- log10(values+1)
  return(res)
}
BoxCox_transformation <- function(values,optimize=TRUE,lambda=c(-5,5),method="Log-Likelihood"){
  require(EnvStats)
  values <- sapply(values,function(x) max(x,1e-10))
  if(optimize){
    boxcox_optimisation <- EnvStats::boxcox(values,lambda=lambda,optimize=optimize,objective.name=method)
    transformed_data <- EnvStats::boxcoxTransform(x=values,lambda=boxcox_optimisation$lambda)
    return(list(data=transformed_data,lambda=boxcox_optimisation$lambda))
  }else{
    transformed_data <- EnvStats::boxcoxTransform(x=values,lambda=lambda)
    return(list(data=transformed_data,lambda=lambda))
  }
}
YeoJohnson_transformation <- function(values,standard=TRUE){
  require(bestNormalize)
  yeo_optimisation <- bestNormalize::yeojohnson(values,standardize=standard)
  return(list(data=yeo_optimisation$x.t,lambda=yeo_optimisation$lambda))
}
Arcsinh_transformation <- function(values,standard=TRUE){
  require(bestNormalize)
  arcsinh_transformation <- bestNormalize::arcsinh_x(values,standardize=standard)
  return(list(data=arcsinh_transformation$x.t))
}
Log_transformation <- function(values,standard=TRUE){
  require(bestNormalize)
  log_transformation <- bestNormalize::log_x(values,standardize=standard,b=10)
  return(list(data=log_transformation$x.t,a=log_transformation$a,b=log_transformation$b))
}

# Data management
Linear_Interpolation_Parameters <- function(Data,markers){
  
  Interpolation_parameters <- NULL
  
  for(m in 1:length(markers)){
    data_markers <- Data[,c("Time","Group","AnimalID",markers[m])]
    colnames(data_markers) <- c("Time","Group","AnimalID","Marker")
    data_markers <- data_markers[!duplicated(data_markers),]
    AnimalIDs <- unique(data_markers$AnimalID)
    for(id in 1:length(AnimalIDs)){
      marker_id <- data_markers[which(data_markers$AnimalID == AnimalIDs[id]),]
      observed_times <- unique(marker_id$Time)
      # Estimation of the parameters (slope and intercept at each observed time points)
      for(t in 1:length(observed_times)){
        t0 <- observed_times[t]
        if(t < length(observed_times)){ # we are not looking at the last observed time point 
          t1 <- observed_times[which(observed_times > t0)][1] # next observed time point
          Y0 <- marker_id$Marker[which(marker_id$Time == t0)]
          Y1 <- marker_id$Marker[which(marker_id$Time == t1)]
          Slope <- (Y0-Y1)/(t0-t1) ; intercept <- Y0 - Slope*t0
        }else{
          Y0 <- marker_id$Marker[which(marker_id$Time == t0)]
          Slope <- 0 ; intercept <- Y0
        }
        tmp_parameters <- data.frame(Time=t0,Group=unique(marker_id$Group),AnimalID=unique(marker_id$AnimalID),Marker=markers[m],
                                     Slope=Slope,Intercept=intercept,stringsAsFactors = FALSE)
        Interpolation_parameters <- rbind(Interpolation_parameters,tmp_parameters,stringsAsFactors = FALSE)
      } # End for t
    } # End for id
  } # End for m
  return(Interpolation_parameters)
}
Combinaison_VL_InterpolationParameters <- function(Data_VL,Parameters,markers){
  New_dataset <- Data_VL
  
  for(m in 1:length(markers)){
    marker <- markers[m]
    tmp_data <- Data_VL
    tmp_data$Slope <- NA ; tmp_data$Intercept <- NA
    tmp_data$Slope <- sapply(seq(1,nrow(tmp_data)),function(i,data,params,marker){
      data_i <- data[i,]
      res <- params$Slope[which(params$Time==data_i$Time & params$Group == data_i$Group & 
                                  params$AnimalID == data_i$AnimalID & params$Marker == marker)]
      return(res)
    },data=tmp_data,params=Parameters,marker=marker)
    tmp_data$Intercept <- sapply(seq(1,nrow(tmp_data)),function(i,data,params,marker){
      data_i <- data[i,]
      res <- params$Intercept[which(params$Time==data_i$Time & params$Group == data_i$Group & 
                                      params$AnimalID == data_i$AnimalID & params$Marker == marker)]
      return(res)
    },data=tmp_data,params=Parameters,marker=marker)
    ind_col <- which(colnames(tmp_data) %in% c("Slope","Intercept"))
    colnames(tmp_data)[ind_col] <- paste(marker,c("Slope","Intercept"),sep=".")
    New_dataset <- cbind(New_dataset,tmp_data[,ind_col],stringsAsFactors = FALSE)
  }# End for m
  
  return(New_dataset)
}

# Covariate selection
Stepwise_CovariateSelection_Adjustment_and_DatasetCreation <- function(step,markers,parameters,fixedCov,data){
  tested_adjustment <- expand.grid(tested_param=parameters,Markers=markers,stringsAsFactors = FALSE)
  tested_adjustment$Label <- tested_adjustment$Markers
  # We remove the markers that have already been selected
  for(i in 1:length(fixedCov)){
    if(!is.null(fixedCov[[i]])){
      tmp_fixed_cov <- fixedCov[[i]]
      if(nrow(tested_adjustment[which(tested_adjustment$tested_param == names(fixedCov)[i] & tested_adjustment$Markers == tmp_fixed_cov),])>0){
        tested_adjustment <- tested_adjustment[-which(tested_adjustment$tested_param == names(fixedCov)[i] & tested_adjustment$Markers == tmp_fixed_cov),] 
      }
    }
  }
  tested_adjustment$Model <- NA
  for(p in 1:length(parameters)){
    tested_adjustment$Model[which(tested_adjustment$tested_param == parameters[p])] <- Model_Files[[parameters[p]]]  
  }
  
  # Reshape of data according to results from previous steps
  tmp <- do.call("rbind",(strsplit(colnames(data),split=".",fixed=TRUE)))
  ind_markers <- which(tmp[,1] == tmp[,2])
  ind_markers <- ind_markers[length(ind_markers)] + 1
  if(!is.null(unlist(fixedCov))){
    New_dataSet <- data[,c(1:(ind_markers-1))]
    for(i in 1:length(fixedCov)){
      if(!is.null(fixedCov[[i]])){
        fixed_cov <- unlist(fixedCov[[i]])
        fixed_cov <- fixed_cov[which(!is.na(fixed_cov))]
        col_ind <- NULL
        col_names <- NULL
        for(j in 1:length(fixed_cov)){
          col_ind <- c(col_ind,which(colnames(data) %in% paste(fixed_cov[j],c("Slope","Intercept"),sep=".")))
          col_names <- c(col_names,paste(fixed_cov[j],c("Slope","Intercept"),names(fixedCov[i]),sep="."))
        }
        Added_columns <- data[,col_ind]
        colnames(Added_columns) <- col_names
        New_dataSet <- cbind(New_dataSet,Added_columns)
      } 
    }
    New_dataSet <- cbind(New_dataSet,data[,ind_markers:ncol(data)])
  }else{
    New_dataSet <- data
  }
  return(list(Dataset=New_dataSet,Adjustment=tested_adjustment))
}
Time_varying_covariate_Monolix <- function(data_file,model_file,Project_file,parameters,Headers,params_distribution,params_variability,interpolation=FALSE){
  Data <- read.table(file=data_file,sep=";",dec=".",header=TRUE,stringsAsFactors=FALSE,check.names = FALSE)
  # Reordering of the monolix headers
  Reorder_Headers <- Headers[intersect(colnames(Data),names(Headers))]
  HeaderTypes <- rep("ignore",ncol(Data)) ; names(HeaderTypes) <- colnames(Data)
  HeaderTypes[which(names(HeaderTypes) %in% names(Reorder_Headers))] <- Reorder_Headers
  
  # Step 1: Creation of the project
  # ----- #
  obs_type <- rep("continuous",4) 
  names(obs_type) <- c("y_1","y_2","y_3","y_4")
  observationNames <- c("y_1","y_2","y_3","y_4")
  # Creation of the Monolix project
  newProject(data=list(dataFile=data_file,
                       header=colnames(Data),
                       headerTypes=as.vector(HeaderTypes),
                       observationNames=observationNames,
                       observationTypes = obs_type,
                       mapping=list('1'="y_1",'2'="y_2",'3'="y_3",'4'="y_4")),
             modelFile=model_file)
  saveProject(projectFile = Project_file)
  
  Scenario <- getScenario()
  Scenario$tasks <- c(populationParameterEstimation=T,conditionalDistributionSampling=T,
                      conditionalModeEstimation=T,standardErrorEstimation=T,
                      logLikelihoodEstimation=T,plots=T)
  
  setScenario(Scenario)
  
  
  # Model Specification
  setErrorModel(y1="constant",y2="constant",y3="constant",y4="constant")
  setIndividualParameterVariability(params_variability)
  setIndividualParameterDistribution(params_distribution)
  
  Init_params <- getPopulationParameterInformation()
  for(p in 1:nrow(Init_params)){
    param <- Init_params$name[p]
    selected_param <- strsplit(param,split="_pop",fixed=TRUE)[[1]]
    if(selected_param %in% parameters$name){
      ind <- which(parameters$name == selected_param)
      Init_params$initialValue[p] <- parameters$values[ind]
      Init_params$method[p] <- parameters$method[ind]
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
}

# Analayze of th results
Results_Extraction_Function <- function(projectFolder,parameter,marker,immuno_param=NA){
  
  results <- data.frame(tested_param=parameter,Markers=marker,LL=NA,BICc=NA,BIC=NA,AIC=NA)
  # likelihood criteria 
  file_LL <- paste(projectFolder,"LogLikelihood", "logLikelihood.txt",sep="/")
  if(file.exists(file_LL)){
    likelihood_criteria <- read.csv(file_LL)
    results$LL <- -0.5*likelihood_criteria$importanceSampling[which(likelihood_criteria$criteria == "-2LL")]
    results$BICc <- likelihood_criteria$importanceSampling[which(likelihood_criteria$criteria == "BICc")]
    results$BIC <- likelihood_criteria$importanceSampling[which(likelihood_criteria$criteria == "BIC")]
    results$AIC <- likelihood_criteria$importanceSampling[which(likelihood_criteria$criteria == "AIC")]
  }
  
  # coefficients of regression linked to immuno markers
  if(!is.na(immuno_param)){
    file_coeff <- paste(projectFolder,"populationParameters.txt",sep="/")
    if(file.exists(file_coeff)){
      pop_parameters <- read.csv(file_coeff)
      coeff_values <- pop_parameters$value[which(pop_parameters$parameter %in% paste(immuno_param,"_pop",sep=""))]
      results <- setNames(cbind(results,coeff_values),c(colnames(results),immuno_param))
      if("se_sa" %in% colnames(pop_parameters)){
        coeff_sd <- pop_parameters$se_sa[which(pop_parameters$parameter %in% paste(immuno_param,"_pop",sep=""))]
        results <- setNames(cbind(results,coeff_sd),c(colnames(results),paste(immuno_param,"_sd",sep="")))
      }else{
        results <- setNames(cbind(results,rep(NA,length(immuno_param))),c(colnames(results),paste(immuno_param,"_sd",sep="")))
      }
    }
  }
  return(results)
}
Random_Effects_variance_Extraction_Function <- function(projectFolder,parameter,marker){
  results <- data.frame(tested_param=parameter,Markers=marker)
  # population parameters (random effects)
  file_coeff <- paste(projectFolder,"populationParameters.txt",sep="/")
  if(file.exists(file_coeff)){
    pop_parameters <- read.csv(file_coeff,stringsAsFactors = FALSE)
    ind <- which(substr(pop_parameters$parameter,1,5) == "omega")
    results <- setNames(cbind(results,t(pop_parameters$value[ind])),c(colnames(results),do.call("rbind",strsplit(pop_parameters$parameter[ind],split="omega_",fixed=T))[,2]))
  }
  return(results)
}
# ----- #


Monolix_Folder <- "Monolix_Estimation"
Rcode_Folder <- "Rcode"
Data_Folder <- "Simulated_data"
Mlxtran_Models_Folder <- "Mlxtran_Models"
Initial_data_file <- paste(Data_Folder,"Simulated_data_ViralLoad_ImmunoMarkers.csv",sep="/")

Folder_CovSelection <- "Linear_Time_varying_covariate"
dir.create(paste(Monolix_Folder,Folder_CovSelection,sep="/"))

Data_Transformation <- "NaturalScale"
dir.create(paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,sep="/"))
# Name of the data file after potential data transformation
New_data_file <- paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,paste("Simulated_data_ViralLoad_ImmunoMarkers_",Data_Transformation,".txt",sep=""),sep="/")


# ------------------------- #
# PART I: Data Management ####

# As input, we assume a dataset shaped as follows (minimum requested)
# AnimalID: column indicating unique identifiant
# Time: column of time of measurements
# log10VL: column gathering observations (whether the type of viral load, gRNA or sgRNA, and the compartment, trachea and nasopharynx)
# obs: column given a unique identifiant to each type of viral load: 1) gRNA in Trachea, 2) sgRNA in Trachea, 3) gRNA in Nasopharynx, 4) sgRNA in Nasopharynx
# Group: column indicating the group of subjects
# Weight: column indicating the weight of animals (used as regressor in the model)  
# censored: column indicating whether observations are left-censored (1) or not (0)
# For each marker, either a column with the time-varying data and/or two columns: one for the intercept and the other for the slope values estimated to re-build linear interpolation


# 0. Download of data ####
Data <- read.table(file=Initial_data_file,dec=".",sep=";",header=TRUE,stringsAsFactors = FALSE,check.names = FALSE)
Markers <- c("ECL","Marker2","Marker3")
ind_markers <- 10 # First column with information about markers (all data of markers are assumed to be gathered at the end of the dataset)

if(length(setdiff(Markers,colnames(Data)))!=0){ # we assume that only linear interpolation parameters are provided
  for(m in 1:length(Markers)){
    marker_data <- Data[,which(colnames(Data) == paste(Markers[m],"Slope",sep="."))]*Data$Time + Data[,which(colnames(Data) == paste(Markers[m],"Intercept",sep="."))]
    Data <- cbind(Data,marker_data)
    colnames(Data)[nrow(Data)] <- Markers[m]
  }
}


# 1. Application of the requested data transformation ####
# If the name of the transformation (variable Data_Transformation) is not one of the transformation indiciated in the "Data_transformation_function",
# no transformation will be applied (e.g. NaturalScale)
if(Data_Transformation %in% c("BoxCox","YeoJohnson-Norm","YeoJohnson","Log")){
  tmp_transformed_data <- Data_transformation_function(Data=Data[,which(colnames(Data) %notin% paste(Markers,rep(c("Slope","Intercept"),each=length(Markers)),sep="."))],markers=Markers,Transf_name=Data_Transformation)
  Reshaped_data <- tmp_transformed_data$data
  Marker_Transf_param <- tmp_transformed_data$lambda
}else{
  tmp_transformed_data <- Data_transformation_function(Data=Data[,which(colnames(Data) %notin% paste(Markers,rep(c("Slope","Intercept"),each=length(Markers)),sep="."))],markers=Markers,Transf_name=Data_Transformation)
  Reshaped_data <- tmp_transformed_data
}

# 2. Estimation of the linear interpolation parameters of the transformed data ####
Parameter_LinInt <- Linear_Interpolation_Parameters(Data=Reshaped_data,markers=Markers)

# 3. Combinaison of VL data & linear interpolation parameters ####
Combined_data <- Combinaison_VL_InterpolationParameters(Data_VL=Reshaped_data[,c(1:(ind_markers-1))],Parameters=Parameter_LinInt,markers=Markers)
# We only keep data with positive time values (after the 2nd exposition)
Combined_data <- Combined_data[which(Combined_data$Time>=0),]
Combined_data <- Combined_data[which(!is.na(Combined_data$log10VL)),]
# Record of the new data
write.table(Combined_data,file=New_data_file,row.names=FALSE,dec=".",sep=";")
# ------------------------- #


# ------------------------- #
# PART II: First Step of the covariate selection ####

# 1. Selection step Settings ####
# ---- #

Step <- "Step1"
Step_folder <- "1TimeCov" 
Tested_parameter <- c("beta","delta") # List of parameters on which adjustment for covariates are tested at this step
dir.create(paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,Step_folder,sep="/"))
for(param in Tested_parameter){
  dir.create(paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,Step_folder,paste("param",param,sep="_"),sep="/"))
}
Model_specification <- paste("TimeVaryingCov_ModelSpecificiation",Step,sep="_") # Name of the file in which the specification of the Monolix models are given for this step (initial values, variability, distribution of parameters, Monolix models to use for the adjustment of covariate on each parameter of interest)
source(paste(Rcode_Folder,paste(Model_specification,".R",sep=""),sep="/"))

# 2. Creation of the dataset and adjustment to test for the step ####
# ---- #

Data_Adjustment_Step1 <- Stepwise_CovariateSelection_Adjustment_and_DatasetCreation(step=Step,markers=Markers,parameters=Tested_parameter,fixedCov = Fixed_timeCov,data=Combined_data)
DataSet_Step1 <- Data_Adjustment_Step1$Dataset
Adjustments_Step1 <- Data_Adjustment_Step1$Adjustment
# Save the new dataset
Data_file <- paste(strsplit(New_data_file,split=".txt",fixed=TRUE)[[1]],"_",Step,".txt",sep="")
write.table(DataSet_Step1,file=Data_file,row.names=FALSE,dec=".",sep=";")

# 3. Model estimations (parallel calculation) ####
# ---- #

ModelEst_Folder <- paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,Step_folder,sep="/")

Nb_Cores <- 19
NCores <- min(detectCores()-1,Nb_Cores)
Cluster <- parallel::makeCluster(NCores,type="SOCK")
registerDoSNOW(Cluster)
pb <- txtProgressBar(max=nrow(Adjustments_Step1[1:nrow(Adjustments_Step1),]),style=3)         # We add a progress bar
progress <- function(n) setTxtProgressBar(pb,n)
opts <- list(progress=progress)

clusterExport(Cluster,c("Adjustments_Step1","Mlxtran_Models_Folder","Rcode_Folder","Model_specification","Data_file","ModelEst_Folder","Time_varying_covariate_Monolix"))


Parallel_results <- foreach(n=seq(1,nrow(Adjustments_Step1)),.errorhandling = "remove",.packages = c("doParallel","foreach","doSNOW","lixoftConnectors","reshape2"),.options.snow = opts)%dopar%{
  initializeLixoftConnectors(software = "monolix", path = "C:/ProgramData/Lixoft/MonolixSuite2019R1")
  adjustment <- Adjustments_Step1[n,]
  
  source(paste(Rcode_Folder,paste(Model_specification,".R",sep=""),sep="/"))
  Project_file <- paste(ModelEst_Folder,paste("param",adjustment$tested_param,sep="_"),paste("TimeCov_",adjustment$tested_param,"_",adjustment$Markers,".mlxtran",sep=""),sep="/")
  # Addition of regressors (linear interpolation parameters) defining the dynamics of the tested marker
  Headers <- c(Header_Monolix,"regressor","regressor") # Addition of 2 regressors, for the slope and the intercept of the tested marker
  names(Headers) <- c(names(Header_Monolix),paste(adjustment$Markers,c("Slope","Intercept"),sep="."))
  
  Time_varying_covariate_Monolix(data_file=Data_file,model_file=paste(Mlxtran_Models_Folder,adjustment$Model,sep="/"),Project_file=Project_file,
                                 parameters=Initial_parameters,Headers=Headers,params_distribution=parameter_distribution,params_variability=parameter_variability)
  
}
close(pb)
stopCluster(Cluster)


# 4. Extraction of the results ####
# ---- #
Step1_projectNames <- paste(ModelEst_Folder,paste("param",Adjustments_Step1$tested_param,sep="_"),paste("TimeCov_",Adjustments_Step1$tested_param,"_",Adjustments_Step1$Markers,sep=""),sep="/")
Step1_Results <- do.call("rbind",lapply(seq(1,length(Step1_projectNames)),function(i) Results_Extraction_Function(Step1_projectNames[i],parameter=Adjustments_Step1$tested_param[i],marker=Adjustments_Step1$Markers[i],immuno_param="beta_immuno")))
Step1_RandomEffect_Sd <- do.call("rbind",lapply(seq(1,length(Step1_projectNames)),function(i) Random_Effects_variance_Extraction_Function(Step1_projectNames[i],parameter=Adjustments_Step1$tested_param[i],marker=Adjustments_Step1$Markers[i])))
# Addition of the results for the models adjusted for group effects or without adjustment
# We extract the results of the models estimated in the R file "Rcode/ODE_model_Estimation_Monolix.R"
Additional_projectNames <- c("Model_Estimation_NoCov","Model_Estimation_CovGroupBeta","Model_Estimation_CovGroupDelta","Model_Estimation_CovGroupBetaDelta")
GroupCov_parameters <- c("None","beta","delta","both")
Additional_project_Results <- do.call("rbind",lapply(seq(1,length(Additional_projectNames)),function(i) Results_Extraction_Function(paste(Monolix_Folder,Additional_projectNames[i],sep="/"),parameter=GroupCov_parameters[i],marker="Group",immuno_param=NA)))
Additional_project_RE_Sd <- do.call("rbind",lapply(seq(1,length(Additional_projectNames)),function(i) Random_Effects_variance_Extraction_Function(paste(Monolix_Folder,Additional_projectNames[i],sep="/"),parameter=GroupCov_parameters[i],marker="Group")))

Global_Step1_Results <- rbind(Step1_Results,cbind(Additional_project_Results,beta_immuno=NA,beta_immuno_sd=NA))
Global_Step1_RE_Sd <- rbind(Step1_RandomEffect_Sd,Additional_project_RE_Sd)
# Record of the results
write.csv(Global_Step1_Results,file=paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,"Step1_criteria_results.csv",sep="/"),row.names = FALSE)
write.csv(Global_Step1_RE_Sd,file=paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,"Step1_RandomEffects_variability.csv",sep="/"),row.names = FALSE)


# 5. Plot of the results ####
# ---- #
require("ggrepel")

Results_Step1 <- read.csv(file=paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,"Step1_criteria_results.csv",sep="/"),header=TRUE,stringsAsFactors = FALSE)
Results_Step1$LL <- -2*Results_Step1$LL
colnames(Results_Step1)[which(colnames(Results_Step1) == "LL")] <- "-2LL"
Results_Step1$RSE <- Results_Step1$beta_immuno_sd*100/abs(Results_Step1$beta_immuno)
Results_Step1$RSE_range <- 1*(Results_Step1$RSE <= 50) + 2*between(Results_Step1$RSE,50,75) + 3*between(Results_Step1$RSE,75,100) + 4*(Results_Step1$RSE>100)
Results_Step1$Num_Marker <- as.numeric(factor(Results_Step1$Marker,levels = unique(Results_Step1$Marker)))
Results_Step1$Color <- 1*(Results_Step1$beta_immuno <= 0) 
Results_Step1$Label <- Results_Step1$Markers

Plot_Step1_BICc <-  ggplot(Results_Step1[which(Results_Step1$Markers != "Group"),]) +
  geom_point(aes(x=Markers,y=BICc,group=Markers,color=as.factor(Color),shape=as.factor(tested_param),size=as.factor(RSE_range)),alpha=0.5,stroke=1.2) + 
  scale_shape_manual(name="Parameter",breaks = c("beta","delta","P","c"),values=c(21,24,22,8),labels=c("beta","delta","P","c")) +
  scale_color_manual(name=paste("Regression","\n","Coefficient",sep=""),breaks=c(0,1),values=c("red","darkblue"),labels=c(">0","<=0")) +
  scale_size_manual(name="RSE (%)", breaks=seq(1,4),values=c(1.5,2,3,5),labels=c("[0,50]","]50,75]","]75,100",">100")) +
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(color="black",size=1.2),
        axis.text.y =  element_text(color="black",size=9),
        axis.text.x =  element_text(color="black",size=9,angle=90,vjust=0.5,hjust=0,face="bold"),
        panel.grid.major.x = element_line(colour="gray90"),
        panel.grid.minor.x = element_blank(),
        axis.title = element_text(color="black",size=12,face="bold"),
        legend.background = element_rect(color="white",size=0.8),
        legend.title = element_text(color="black",size=10,face="bold"),
        legend.text = element_text(color="black",face="italic",size=10),
        legend.key.width = unit(1.2,"cm"), legend.key.height = unit(0.5,"cm"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        strip.background = element_rect(fill="white",color="black"),
        strip.text = element_text(color="black",size=10),
        panel.spacing.x=unit(0, "lines"),
        axis.ticks.x = element_line(color="black"),
        panel.border = element_rect(color="black",fill=NA,size=1.2),
        legend.position = "right",
        legend.direction = "vertical",
        legend.justification = "center",
        legend.box = "vertical") +
  scale_x_discrete(position="top") + 
  scale_y_continuous(breaks = c(seq(500,1000,by=2))) +
  xlab("") + ylab("BICc") +
  guides(color=guide_legend(override.aes = list(size=3)),shape=guide_legend(override.aes = list(size=3))) + 
  # Addition of the lines of reference
  geom_hline(aes(yintercept=Results_Step1$BICc[which(Results_Step1$Markers == "Group" & Results_Step1$tested_param == "None")],linetype="No_Cov"),color="gray40",size=1,alpha=0.75) + 
  geom_hline(aes(yintercept=Results_Step1$BICc[which(Results_Step1$Markers == "Group" & Results_Step1$tested_param == "beta")],linetype="Group_CovBeta"),color="darkblue",size=1,alpha=0.75) + 
  geom_hline(aes(yintercept=Results_Step1$BICc[which(Results_Step1$Markers == "Group" & Results_Step1$tested_param == "delta")],linetype="Group_CovDelta"),color="forestgreen",size=1,alpha=0.75) + 
  geom_hline(aes(yintercept=Results_Step1$BICc[which(Results_Step1$Markers == "Group" & Results_Step1$tested_param == "both")],linetype="Group_CovBetaDelta"),color="darkorange",size=1,alpha=0.75) + 
  # Addition of Labels
  geom_label_repel(aes(x=Markers,y=BICc,label=Label),force=5,color="darkorange",size=3,alpha=1) +
  scale_linetype_manual(name="References",breaks = c("No_Cov","Group_CovBeta","Group_CovDelta","Group_CovBetaDelta"),values=c("solid","longdash","twodash","dotted"),
                        labels=c(paste("Model without covariates",sep=""),
                                 paste("Model with 'Group'","\n","covariate on beta",sep=""),
                                 paste("Model with 'Group'","\n","covariate on delta",sep=""),
                                 paste("Model with 'Group'","\n","covariates on beta & delta",sep="")),drop=FALSE) +
  guides(linetype = guide_legend(ncol=2,override.aes = list(color = c("gray40","darkblue","forestgreen","darkorange"))),
         size = guide_legend(ncol=2)) +
  theme(legend.position = "bottom",legend.direction = "vertical",legend.box = "horizontal") 
  

# 6. Replacement of the group Effect ? ####
# For the model displaying the lowest value of BICc (here), we adjust the model on the group effect, in addition to the marker
Best_model_Step1 <- Results_Step1[which(Results_Step1$BICc == min(Results_Step1$BICc[which(Results_Step1$Markers != "Group")],na.rm=TRUE)),]
# As excepted by the simulation, we select the marker ECL on the parameter beta
# We load the mononolix project
Best_model_Step1_name <- paste(ModelEst_Folder,paste("param",Best_model_Step1$tested_param,sep="_"),paste("TimeCov_",Best_model_Step1$tested_param,"_",Best_model_Step1$Markers,sep=""),sep="/")
loadProject(projectFile = paste(Best_model_Step1_name,".mlxtran",sep=""))
saveProject(projectFile = paste(Best_model_Step1_name,paste("_CovGroup",Best_model_Step1$tested_param,sep=""),".mlxtran",sep=""))
# Addition of the group covariates
addCategoricalTransformedCovariate(groupVsnaive=list(reference="G_Naive",from="Group",transformed=list(G_Naive="Naive",G_conv="Convalescent",G_convCD40="Conv-CD40")))
ifelse(Best_model_Step1$tested_param == "beta",setCovariateModel(beta_N_pow=c(groupVsnaive=TRUE)),setCovariateModel(delta_N=c(groupVsnaive=TRUE)))
saveProject()
# Model Estimation
runScenario()
saveProject()


# 7. Explained variability ####
RandomEffect_Variability_Step1 <- read.csv(file=paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,"Step1_RandomEffects_variability.csv",sep="/"),header=TRUE,stringsAsFactors = FALSE)
# The model without covariates is used as reference for the estimation of explained variability
RandomEffect_Variability_Step1$Expl.Var.beta <- 100 - RandomEffect_Variability_Step1$beta_N_pow*100/RandomEffect_Variability_Step1$beta_N_pow[which(RandomEffect_Variability_Step1$tested_param == "None")]
RandomEffect_Variability_Step1$Expl.Var.delta <- 100 - RandomEffect_Variability_Step1$delta_N*100/RandomEffect_Variability_Step1$delta_N[which(RandomEffect_Variability_Step1$tested_param == "None")]


# PART III: Second Step of the covariate selection ####
# 1. Selection step Settings ####
Step <- "Step2"
Step_folder <- "2TimeCov_betaECLRBD" 
Tested_parameter <- c("beta","delta") # List of parameters on which adjustment for covariates are tested at this step
dir.create(paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,Step_folder,sep="/"))
for(param in Tested_parameter){
  dir.create(paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,Step_folder,paste("param",param,sep="_"),sep="/"))
}
Combined_data <- read.table(file=New_data_file,header=TRUE,dec=".",sep=";")
Model_specification <- paste("TimeVaryingCov_ModelSpecificiation",Step,sep="_") # Name of the file in which the specification of the Monolix models are given for this step (initial values, variability, distribution of parameters, Monolix models to use for the adjustment of covariate on each parameter of interest)
source(paste(Rcode_Folder,paste(Model_specification,".R",sep=""),sep="/"))


# 2. Creation of the dataset and adjustment to test for the step ####
# ---- #
Combined_data <- read.table(file=New_data_file,header=TRUE,dec=".",sep=";")

Data_Adjustment_Step2 <- Stepwise_CovariateSelection_Adjustment_and_DatasetCreation(step=Step,markers=Markers,parameters=Tested_parameter,fixedCov = Fixed_timeCov,data=Combined_data)
DataSet_Step2 <- Data_Adjustment_Step2$Dataset
Adjustments_Step2 <- Data_Adjustment_Step2$Adjustment
# Modification of the Monolix headers
if(is.null(unlist(Fixed_timeCov))){
  Header_Monolix <- c("time","catcov","id","regressor","obsid","observation","cens")
  names(Header_Monolix) <- c("Time","Group","AnimalID","Weight","obs","log10VL","censored")
}else{
  Header_Monolix <- c("time","catcov","id","regressor","obsid","observation","cens",rep("regressor",2*length(unlist(Fixed_timeCov))))
  names(Header_Monolix) <- c("Time","Group","AnimalID","Weight","obs","log10VL","censored",colnames(DataSet_Step2)[10:(9+2*length(unlist(Fixed_timeCov)))])
}
# Save the new dataset
Data_file <- paste(strsplit(New_data_file,split=".txt",fixed=TRUE)[[1]],"_",Step,".txt",sep="")
write.table(DataSet_Step2,file=Data_file,row.names=FALSE,dec=".",sep=";")


# 3. Model estimations (parallel calculation) ####
# ---- #

ModelEst_Folder <- paste(Monolix_Folder,Folder_CovSelection,Data_Transformation,Step_folder,sep="/")

Nb_Cores <- 5
NCores <- min(detectCores()-1,Nb_Cores)
Cluster <- parallel::makeCluster(NCores,type="SOCK")
registerDoSNOW(Cluster)
pb <- txtProgressBar(max=nrow(Adjustments_Step2[1:nrow(Adjustments_Step2),]),style=3)         # We add a progress bar
progress <- function(n) setTxtProgressBar(pb,n)
opts <- list(progress=progress)

clusterExport(Cluster,c("Adjustments_Step2","Mlxtran_Models_Folder","Rcode_Folder","Model_specification","Header_Monolix","Data_file","ModelEst_Folder","Time_varying_covariate_Monolix"))

# n <- 1
Parallel_results <- foreach(n=seq(1,nrow(Adjustments_Step2)),.errorhandling = "remove",.packages = c("doParallel","foreach","doSNOW","lixoftConnectors","reshape2"),.options.snow = opts)%dopar%{
  initializeLixoftConnectors(software = "monolix", path = "C:/ProgramData/Lixoft/MonolixSuite2019R1")
  adjustment <- Adjustments_Step2[n,]
  
  source(paste(Rcode_Folder,paste(Model_specification,".R",sep=""),sep="/"))
  Project_file <- paste(ModelEst_Folder,paste("param",adjustment$tested_param,sep="_"),paste("TimeCov_",adjustment$tested_param,"_",adjustment$Markers,".mlxtran",sep=""),sep="/")
  # Addition of regressors (linear interpolation parameters) defining the dynamics of the tested marker
  Headers <- c(Header_Monolix,"regressor","regressor") # Addition of 2 regressors, for the slope and the intercept of the tested marker
  names(Headers) <- c(names(Header_Monolix),paste(adjustment$Markers,c("Slope","Intercept"),sep="."))
  
  Time_varying_covariate_Monolix(data_file=Data_file,model_file=paste(Mlxtran_Models_Folder,adjustment$Model,sep="/"),Project_file=Project_file,
                                 parameters=Initial_parameters,Headers=Headers,params_distribution=parameter_distribution,params_variability=parameter_variability)
  
}
close(pb)
stopCluster(Cluster)
