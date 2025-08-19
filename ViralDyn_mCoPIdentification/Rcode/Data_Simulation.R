# ----------- #
# OBJECTIVE: Simulation of viral load and immuno markers dynamics 
# Author: Marie Alexandre
# Date: 2021/11/23

# R version: 3.6.2
# Simulx version: 2019R1

# --- Description: 
# In this file, we aim at simulated data that will be used by the other files and thus providing an overview of the use of the codes without raw data. 
# Viral loads (gRNA and sgRNA) are simulated via the mechanistic model based on ODE with immuno markers considered as external variable (covariate)
# Viral load within trachea and nasopharynx are jointly simulated via the mlxR library while immuno markers mimicking the shape of pseudo-neutralization data (see Figure 4A) are simulated via logistic functions

# Tips: Ctrl + shift + O to extend document outline
# ----------- #


rm(list=ls())
'%notin%' <- Negate('%in%')


# Libraries --- ####
require(mlxR) ; require(lixoftConnectors)
require(plyr)
# ----- #

# Initialization of the lixoft connectors (for simulx) --- ####
Monolix_path <- "C:/ProgramData/Lixoft/MonolixSuite2019R1"  # To modify if necessary
initializeLixoftConnectors(software = "simulx", path = Monolix_path,force=TRUE)
initMlxR(path=Monolix_path)
# ----- #

# Functions --- ####
# 1. Functions for the simulation of immuno markers ####
Linear_FUNCTION <- function(t,a,b){
  return(a*t+b)
}
Logistic_FUNCTION <- function(t,R0,Rf,k,x0){
  return((R0-Rf)/(1+exp(-k*(x0-t)))+Rf)
}
Simulation_log10ECL_WithoutNoise_FUNCTION <- function(time,parameters,Group,Nb_subject){
  
  simulated_ECL <- expand.grid(Time=time,AnimalID=seq(1,Nb_subject,by=1),Group=Group,log10ECL=NA)
  
  if(Group == "Convalescent"){ # Logistic function for the log10 transformation of ECL
    simulated_parameters <- expand.grid(AnimalID=seq(1,Nb_subject,by=1),Group=Group,R0=NA,Rf=NA,k=NA,x0=NA)
    # simulation of the random effect on the slope (a) and the intercept (b)
    eta_R0 <- rnorm(Nb_subject,mean=0,sd=parameters["sd_R0"])
    eta_Rf <- rnorm(Nb_subject,mean=0,sd=parameters["sd_Rf"])
    eta_k <- rnorm(Nb_subject,mean=0,sd=parameters["sd_k"])
    simulated_parameters$R0 <- parameters["R0"] + eta_R0
    simulated_parameters$Rf <- parameters["Rf"] + eta_Rf
    simulated_parameters$k <- parameters["k"] + eta_k
    simulated_parameters$x0 <- parameters["x0"] 
    
    simulated_ECL$log10ECL <- sapply(seq(1,nrow(simulated_ECL)),function(i,params){
      res <- Logistic_FUNCTION(t=simulated_ECL$Time[i],R0=params$R0[which(params$AnimalID == simulated_ECL$AnimalID[i])],
                               Rf=params$Rf[which(params$AnimalID == simulated_ECL$AnimalID[i])],
                               k=params$k[which(params$AnimalID == simulated_ECL$AnimalID[i])],
                               x0=params$x0[which(params$AnimalID == simulated_ECL$AnimalID[i])])
      return(res)
    },params=simulated_parameters)
  }else{ # Linear function for the log10 transformation of ECL
    simulated_parameters <- expand.grid(AnimalID=seq(1,Nb_subject,by=1),Group=Group,a=NA,b=NA)
    # simulation of the random effect on the slope (a) and the intercept (b)
    eta_a <- rnorm(Nb_subject,mean=0,sd=parameters["sd_a"])
    eta_b <- rnorm(Nb_subject,mean=0,sd=parameters["sd_b"])
    simulated_parameters$a <- parameters["a"] + eta_a
    simulated_parameters$b <- parameters["b"] + eta_b
    
    simulated_ECL$log10ECL <- sapply(seq(1,nrow(simulated_ECL)),function(i,params){
      res <- Linear_FUNCTION(t=simulated_ECL$Time[i],a=params$a[which(params$AnimalID == simulated_ECL$AnimalID[i])],b=params$b[which(params$AnimalID == simulated_ECL$AnimalID[i])])
    },params=simulated_parameters)
  }
  
  return(list(Simulated_ECL=simulated_ECL,Simulated_params=simulated_parameters))
}

# 2. Functions for the simulation of viral load ####
Linear_Interpolation_Parameters <- function(Data,markers){
  # Function estimating the parameters for the linear interpolation of the data (slope and intercept at each observed time point)
  
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
Simulation_VL_FUNCTION <- function(parameters,time_simulation,ecl_data,model){
  # In this function, the dataset ecl_data must contains the linear interpolation parameters at each time point
  
  Simulated_ViralLoad <- NULL
  Simulated_parameters <- NULL
  animalIDs <- unique(ecl_data$AnimalID)
  for(i in 1:length(animalIDs)){
    # i <- 1
    ecl_data_id <- ecl_data[which(ecl_data$AnimalID == animalIDs[i]),]
    weight_id <- unique(ecl_data_id$Weight)
    group_i <- unique(ecl_data_id$Group)
    
    ## SIMULATION 
    # Outputs
    f <- list(name=c("f_VLTrachea","f_VLSGTrachea","f_VLNaso","f_VLSGNaso"),time=time_simulation)
    y <- list(name=c("VLTrachea","VLSGTrachea","VLNaso","VLSGNaso"),time=time_simulation)
    M <- list(name=c("Marker"),time=time_simulation)
    g <- list(size=1,level='individual')
    
    # Regressors
    reg_weight <- list(name="Weight",time=time_simulation,value=weight_id)
    reg_slope <- list(name="Marker_slope",time=unique(ecl_data_id$Time),value=ecl_data_id$ECL.Slope)
    reg_int <- list(name="Marker_intercept",time=unique(ecl_data_id$Time),value=ecl_data_id$ECL.Intercept)
    regressors <- list(reg_weight,reg_slope,reg_int)
    
    simulation_results <- simulx(model=model,parameter = parameters,output = list(f,y,M,list(name=c("beta_N_pow","delta_N"))),regressor=regressors,group=g)
    
    # Addition of the linear interpolation parameters
    linear_interpolation_params <- Linear_Interpolation_Parameters(Data=data.frame(Time=time_simulation,Group=group_i,AnimalID=animalIDs[i],Marker=simulation_results$Marker$Marker),markers = "Marker")
    
    # Simulated viral load
    VLTrachea_noise <- data.frame(AnimalID=animalIDs[i],Time=time_simulation,
                                  log10VL=simulation_results$VLTrachea$VLTrachea,obs=1,Group=group_i,type="Total",
                                  organ="Trachea",Weight=weight_id,ECL=simulation_results$Marker$Marker,
                                  ECL.Slope=linear_interpolation_params$Slope,ECL.Intercept=linear_interpolation_params$Intercept)
    VLSGTrachea_noise <- data.frame(AnimalID=animalIDs[i],Time=time_simulation,
                                    log10VL=simulation_results$VLSGTrachea$VLSGTrachea,obs=2,Group=group_i,type="Sg",
                                    organ="Trachea",Weight=weight_id,ECL=simulation_results$Marker$Marker,
                                    ECL.Slope=linear_interpolation_params$Slope,ECL.Intercept=linear_interpolation_params$Intercept)
    VLNaso_noise <- data.frame(AnimalID=animalIDs[i],Time=time_simulation,
                               log10VL=simulation_results$VLNaso$VLNaso,obs=3,Group=group_i,type="Total",
                               organ="Naso",Weight=weight_id,ECL=simulation_results$Marker$Marker,
                               ECL.Slope=linear_interpolation_params$Slope,ECL.Intercept=linear_interpolation_params$Intercept)
    VLSGNaso_noise <- data.frame(AnimalID=animalIDs[i],Time=time_simulation,
                                 log10VL=simulation_results$VLSGNaso$VLSGNaso,obs=4,Group=group_i,type="Sg",
                                 organ="Naso",Weight=weight_id,ECL=simulation_results$Marker$Marker,
                                 ECL.Slope=linear_interpolation_params$Slope,ECL.Intercept=linear_interpolation_params$Intercept)
    
    simulated_viralload_id <- rbind(VLTrachea_noise,VLSGTrachea_noise,VLNaso_noise,VLSGNaso_noise,stringsAsFactors=FALSE)
    Simulated_ViralLoad <- rbind(Simulated_ViralLoad,simulated_viralload_id,stringsAsFactors=FALSE)
    
    Simulated_parameters <- rbind(Simulated_parameters,data.frame(AnimalID=animalIDs[i],Group=group_i,simulation_results$parameter))
  }
  return(list(ViralLoad=Simulated_ViralLoad,Parameters=Simulated_parameters))
}
# ----- #

# Model of simulation --- ####
Model_file <- "Rcode/MlxtranModel_MechanisticModel_SARSCoV2_VL_Simulation.txt"
# ----- #

# ------------------------- #
# PART 0: Definition of features of simulation ####
# -- Timepoints of simulation (based on the main study [Marlin, 2021])
Time_VL <- c(0,1,2,3,4,6,9,14,20)   # Time of final observations (gRNA)
Time_sgVL <- c(0,1,2,3,4,6)         # Time of final observations (sgRNA)
Time_markers <- c(0,4,9,20)         # Time of final observations (markers, pseudo-neutralization)
Time_simulation <- seq(0,20,by=1)   # Time used for simulation
# -- Information about groups and animals 
Nb_animal_group <- 6 # Number of animals by group
Group_Names <- c("Naive","Convalescent","Conv-CD40")  # WARNING: the name of the group is used in the function of marker simulation to determine the shape of the dynamics
# -- Left-censorship (based on the main study [Marlin, 2021])
Cens_VL <- log10(476)    # left-censoring threshold for gRNA load (natural scale)
Cens_sgVL <- log10(749)  # left-censoring threshold for sgRNA load (natural scale)

Nb_markers <- 3 # number of immuno markers simulated to mimick the pool of markers in raw data (only one marker will be used as external covariate in the simulation of viral load to be more linked with viral load than others)
SD_noise_markers <- c(0.2,0.3,0.3) # standard deviation of the white noise added to simulated markers to obtain unperfect data

Simulation_folder <- "Simulated_data"
# ------------------------- #




# ------------------------- #
# PART I: Simulation of immuno markers (based on the shape of pseudo-neutralization = logistic and linear functions according to groups) ####

# > 0. Simulation of animal weight (uniform distribution based on raw data)
Simulated_weights <- data.frame(AnimalID=c(rep(seq(1,Nb_animal_group),3)+rep(c(100,200,300),each=Nb_animal_group)),
                                Group=rep(Group_Names,each=Nb_animal_group),
                                Weight=runif(n=length(Group_Names)*Nb_animal_group,min=3.00,max=6.230))

# > 1. Simulation of the marker the optimal marker selected on raw data: binding inhibition between RBD domain and ACE2 receptor (labelled ECL) ####

# > 1.a. Simulation of marker dynamics without noise
# Individual dynamics of immuno markers are firstly simulated without noise (to be used in the simulation of the viral load)
# Model parameters used to simulate the data were estimated on raw data
Simulated_log10ECL_Naive <- Simulation_log10ECL_WithoutNoise_FUNCTION(time=Time_markers,parameters=c(a=-0.02,b=5.557,sd_a=0.01,sd_b=0.15),Group="Naive",Nb_subject=Nb_animal_group)
Simulated_log10ECL_Conv <- Simulation_log10ECL_WithoutNoise_FUNCTION(time=Time_markers,parameters=c(R0=5.2,Rf=1.7,k=1.5,x0=7,sd_R0=0.048,sd_Rf=0.09,sd_k=0.075),Group="Convalescent",Nb_subject=Nb_animal_group)
Simulated_log10ECL_Vacc <- Simulation_log10ECL_WithoutNoise_FUNCTION(time=Time_markers,parameters=c(a=-0.004,b=1.75,sd_a=0.002,sd_b=0.1),Group="Conv-CD40",Nb_subject=Nb_animal_group)
# Combination of the simulated data
Simulated_log10ECL <- rbind(Simulated_log10ECL_Naive$Simulated_ECL,Simulated_log10ECL_Conv$Simulated_ECL,Simulated_log10ECL_Vacc$Simulated_ECL)
Simulated_log10ECL$AnimalID <- Simulated_log10ECL$AnimalID + 100*(Simulated_log10ECL$Group == "Naive") + 200*(Simulated_log10ECL$Group == "Convalescent") + 300*(Simulated_log10ECL$Group == "Conv-CD40")
Simulated_parameters <- plyr::rbind.fill(Simulated_log10ECL_Naive$Simulated_params,Simulated_log10ECL_Conv$Simulated_params,Simulated_log10ECL_Vacc$Simulated_params)
Simulated_parameters$AnimalID <- Simulated_parameters$AnimalID + 100*(Simulated_parameters$Group == "Naive") + 200*(Simulated_parameters$Group == "Convalescent") + 300*(Simulated_parameters$Group == "Conv-CD40")
# Addition of Weight in the dataset
Simulated_log10ECL$Weight <- sapply(seq(1,nrow(Simulated_log10ECL)),function(i) return(Simulated_weights$Weight[which(Simulated_weights$AnimalID == Simulated_log10ECL$AnimalID[i])]))
Simulated_log10ECL <- Simulated_log10ECL[,c("Time","AnimalID","Group","Weight","log10ECL")]

# > 1.b. Simulation of marker dynamics with noise
# The simulated white noise is added on log10-transformed dynamics
Simulated_noise <- rnorm(nrow(Simulated_log10ECL),mean=0,sd=SD_noise_markers[1])
Simulated_ECL_Noised <- Simulated_log10ECL 
Simulated_ECL_Noised$log10ECL <- Simulated_ECL_Noised$log10ECL + Simulated_noise
Simulated_ECL_Noised$ECL <- 10^(Simulated_ECL_Noised$log10ECL)

Simulated_log10ECL$ECL <- 10^(Simulated_log10ECL$log10ECL)
# Record of the individual parameters simulated for the first marker and the simulated data
write.table(Simulated_log10ECL,file=paste(Simulation_folder,paste("Simulated_ECLRBD_ExternalData_UnNoised.csv",sep=""),sep="/"),row.names=FALSE,sep=";")
write.table(Simulated_ECL_Noised,file=paste(Simulation_folder,paste("Simulated_ECLRBD_ExternalData_Noised.csv",sep=""),sep="/"),row.names=FALSE,sep=";")
write.table(Simulated_parameters,file=paste(Simulation_folder,paste("Simulated_parameters_ECLRBD_ExternalData.csv",sep=""),sep="/"),row.names=FALSE,sep=";")


# > 2. Simulation of the other markers  ####
# As 2nd marker, we simulate a marker with the same dynamics between convalescent and conv-vaccinated animals (logistic) and we assume a larger variance of white noise
Simulated_log10Marker2_Naive <- Simulation_log10ECL_WithoutNoise_FUNCTION(time=Time_markers,parameters=c(a=-0.02,b=5.557,sd_a=0.01,sd_b=0.15),Group="Naive",Nb_subject=Nb_animal_group)
Simulated_log10Marker2_Conv <- Simulation_log10ECL_WithoutNoise_FUNCTION(time=Time_markers,parameters=c(R0=5.2,Rf=1.7,k=1.5,x0=7,sd_R0=0.048,sd_Rf=0.09,sd_k=0.075),Group="Convalescent",Nb_subject=Nb_animal_group)
Simulated_log10Marker2_Vacc <- Simulation_log10ECL_WithoutNoise_FUNCTION(time=Time_markers,parameters=c(R0=5.2,Rf=1.7,k=1.5,x0=7,sd_R0=0.048,sd_Rf=0.09,sd_k=0.075),Group="Convalescent",Nb_subject=Nb_animal_group)
# We modify the name of the group for vaccinated animals
Simulated_log10Marker2_Vacc$Simulated_ECL$Group <- "Conv-CD40"
Simulated_log10Marker2_Vacc$Simulated_params$Group <- "Conv-CD40"
# Combination of the simulated data
Simulated_log10Marker2 <- rbind(Simulated_log10Marker2_Naive$Simulated_ECL,Simulated_log10Marker2_Conv$Simulated_ECL,Simulated_log10Marker2_Vacc$Simulated_ECL)
Simulated_log10Marker2$AnimalID <- Simulated_log10Marker2$AnimalID + 100*(Simulated_log10Marker2$Group == "Naive") + 200*(Simulated_log10Marker2$Group == "Convalescent") + 300*(Simulated_log10Marker2$Group == "Conv-CD40")
Simulated_parameters_Marker2 <- plyr::rbind.fill(Simulated_log10Marker2_Naive$Simulated_params,Simulated_log10Marker2_Conv$Simulated_params,Simulated_log10Marker2_Vacc$Simulated_params)
Simulated_parameters_Marker2$AnimalID <- Simulated_parameters_Marker2$AnimalID + 100*(Simulated_parameters_Marker2$Group == "Naive") + 200*(Simulated_parameters_Marker2$Group == "Convalescent") + 300*(Simulated_parameters_Marker2$Group == "Conv-CD40")
# Addition of noise
Simulated_noise_Marker2 <- rnorm(nrow(Simulated_log10Marker2),mean=0,sd=SD_noise_markers[2])
Simulated_log10Marker2_Noised <- Simulated_log10Marker2 
Simulated_log10Marker2_Noised$log10ECL <- Simulated_log10Marker2_Noised$log10ECL + Simulated_noise_Marker2
Simulated_log10Marker2_Noised$ECL <- 10^(Simulated_log10Marker2_Noised$log10ECL)
Simulated_log10Marker2$ECL <- 10^(Simulated_log10Marker2$log10ECL)
# Record of the individual parameters simulated and the simulated data
write.table(Simulated_log10Marker2,file=paste(Simulation_folder,paste("Simulated_Marker2_ExternalData_UnNoised.csv",sep=""),sep="/"),row.names=FALSE,sep=";")
write.table(Simulated_log10Marker2_Noised,file=paste(Simulation_folder,paste("Simulated_Marker2_ExternalData_Noised.csv",sep=""),sep="/"),row.names=FALSE,sep=";")
write.table(Simulated_parameters_Marker2,file=paste(Simulation_folder,paste("Simulated_parameters_Marker2_ExternalData.csv",sep=""),sep="/"),row.names=FALSE,sep=";")


# As 3rd marker, we simulate a marker with the same dynamics (linear) between the 3 groups
Simulated_log10Marker3_Naive <- Simulation_log10ECL_WithoutNoise_FUNCTION(time=Time_markers,parameters=c(a=-0.02,b=5.557,sd_a=0.01,sd_b=0.15),Group="Naive",Nb_subject=Nb_animal_group)
Simulated_log10Marker3_Conv <- Simulation_log10ECL_WithoutNoise_FUNCTION(time=Time_markers,parameters=c(a=-0.02,b=5.557,sd_a=0.01,sd_b=0.15),Group="Naive",Nb_subject=Nb_animal_group)
Simulated_log10Marker3_Vacc <- Simulation_log10ECL_WithoutNoise_FUNCTION(time=Time_markers,parameters=c(a=-0.02,b=5.557,sd_a=0.01,sd_b=0.15),Group="Naive",Nb_subject=Nb_animal_group)
# We modify the name of the group for vaccinated and convalescent animals
Simulated_log10Marker3_Conv$Simulated_ECL$Group <- "Convalescent"
Simulated_log10Marker3_Conv$Simulated_params$Group <- "Convalescent"
Simulated_log10Marker3_Vacc$Simulated_ECL$Group <- "Conv-CD40"
Simulated_log10Marker3_Vacc$Simulated_params$Group <- "Conv-CD40"
# Combination of the simulated data
Simulated_log10Marker3 <- rbind(Simulated_log10Marker3_Naive$Simulated_ECL,Simulated_log10Marker3_Conv$Simulated_ECL,Simulated_log10Marker3_Vacc$Simulated_ECL)
Simulated_log10Marker3$AnimalID <- Simulated_log10Marker3$AnimalID + 100*(Simulated_log10Marker3$Group == "Naive") + 200*(Simulated_log10Marker3$Group == "Convalescent") + 300*(Simulated_log10Marker3$Group == "Conv-CD40")
Simulated_parameters_Marker3 <- plyr::rbind.fill(Simulated_log10Marker3_Naive$Simulated_params,Simulated_log10Marker3_Conv$Simulated_params,Simulated_log10Marker3_Vacc$Simulated_params)
Simulated_parameters_Marker3$AnimalID <- Simulated_parameters_Marker3$AnimalID + 100*(Simulated_parameters_Marker3$Group == "Naive") + 200*(Simulated_parameters_Marker3$Group == "Convalescent") + 300*(Simulated_parameters_Marker3$Group == "Conv-CD40")
# Addition of noise
Simulated_noise_Marker3 <- rnorm(nrow(Simulated_log10Marker3),mean=0,sd=SD_noise_markers[2])
Simulated_log10Marker3_Noised <- Simulated_log10Marker3 
Simulated_log10Marker3_Noised$log10ECL <- Simulated_log10Marker3_Noised$log10ECL + Simulated_noise_Marker3
Simulated_log10Marker3_Noised$ECL <- 10^(Simulated_log10Marker3_Noised$log10ECL)
Simulated_log10Marker3$ECL <- 10^(Simulated_log10Marker3$log10ECL)
# Record of the individual parameters simulated and the simulated data
write.table(Simulated_log10Marker3,file=paste(Simulation_folder,paste("Simulated_Marker3_ExternalData_UnNoised.csv",sep=""),sep="/"),row.names=FALSE,sep=";")
write.table(Simulated_log10Marker3_Noised,file=paste(Simulation_folder,paste("Simulated_Marker3_ExternalData_Noised.csv",sep=""),sep="/"),row.names=FALSE,sep=";")
write.table(Simulated_parameters_Marker3,file=paste(Simulation_folder,paste("Simulated_parameters_Marker3_ExternalData.csv",sep=""),sep="/"),row.names=FALSE,sep=";")
# ------------------------- #



# ------------------------- #
# PART II: Simulation of viral load using simulated ECL as external covariate (UnNoised) ####
Simulated_ECL <- read.table(file=paste(Simulation_folder,paste("Simulated_ECLRBD_ExternalData_UnNoised.csv",sep=""),sep="/"),header=TRUE,stringsAsFactors = FALSE,sep=";")

# > 1. Linear Interpolation of the pseudo-neutralization data ####
# In our approach, markers are used as continuous time-varying covariates using a linear interpolation of the data
# Similarly, we need to use linearly interpolated dynamics of markers to simulated viral dynamics
Interpolation_param_ECL <- Linear_Interpolation_Parameters(Data=Simulated_ECL,markers="ECL")
Simulated_ECL$ECL.Slope <- sapply(seq(1,nrow(Simulated_ECL)),function(i,param,data){
  tmp <- param$Slope[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=Interpolation_param_ECL,data=Simulated_ECL)
Simulated_ECL$ECL.Intercept <- sapply(seq(1,nrow(Simulated_ECL)),function(i,param,data){
  tmp <- param$Intercept[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=Interpolation_param_ECL,data=Simulated_ECL)

# > 2. Simulation of Viral load ####
# Definition of Model parameter (based on estimations obtained on real data)
Naive_parameters <- c(beta_N_pow_pop=-8.1,fact_beta_T=0,delta_N_pop=1.0,fact_delta_T=0,c=3,cI=20,k=3,mu=1e-3,g=0,
                      P_N=1e4,fact_P_T=-2.64,alpha_VLSG=1.33,thresh_Weight=4.5,phi_immuno=1.74e-5,
                      omega_beta_N_pow=2.23e-1,omega_delta_N=1.96e-1,a1=1.25,a2=1.33,a3=1.13,a4=1.57)

Convalescent_parameters <- c(beta_N_pow_pop=-8.1,fact_beta_T=0,delta_N_pop=1.0*exp(0.575),fact_delta_T=0,c=3,cI=20,k=3,mu=1e-3,g=0,
                             P_N=1e4,fact_P_T=-2.64,alpha_VLSG=1.33,thresh_Weight=4.5,phi_immuno=1.74e-5,
                             omega_beta_N_pow=2.23e-1,omega_delta_N=1.96e-1,a1=1.25,a2=1.33,a3=1.13,a4=1.57)

Vaccinated_parameters <- c(beta_N_pow_pop=-8.1,fact_beta_T=0,delta_N_pop=1.0*exp(0.778),fact_delta_T=0,c=3,cI=20,k=3,mu=1e-3,g=0,
                           P_N=1e4,fact_P_T=-2.64,alpha_VLSG=1.33,thresh_Weight=4.5,phi_immuno=1.74e-5,
                           omega_beta_N_pow=2.23e-1,omega_delta_N=1.96e-1,a1=1.25,a2=1.33,a3=1.13,a4=1.57)

# Simulation of Viral load (based on ecl without noise)
# Observation id: 1) gRNA in trachea, 2) sgRNA in Trachea, 3) gRNA in Nasopharynx, 4) sgRNA in Nasopharynx
Simulation_Naive <- Simulation_VL_FUNCTION(parameters=Naive_parameters,time_simulation=Time_simulation,ecl_data=Simulated_ECL[which(Simulated_ECL$Group == "Naive"),],model=Model_file)
Simulation_Conv <- Simulation_VL_FUNCTION(parameters=Convalescent_parameters,time_simulation=Time_simulation,ecl_data=Simulated_ECL[which(Simulated_ECL$Group == "Convalescent"),],model=Model_file)
Simulation_Vacc <- Simulation_VL_FUNCTION(parameters=Vaccinated_parameters,time_simulation=Time_simulation,ecl_data=Simulated_ECL[which(Simulated_ECL$Group == "Conv-CD40"),],model=Model_file)

Simulated_ViralLoad <- rbind(Simulation_Naive$ViralLoad,Simulation_Conv$ViralLoad,Simulation_Vacc$ViralLoad)
Simulated_Parameters <- rbind(Simulation_Naive$Parameters,Simulation_Conv$Parameters,Simulation_Vacc$Parameters)
# Addition of censoring information
Simulated_ViralLoad$censored <- 1*(Simulated_ViralLoad$obs %in% c(1,3))*(Simulated_ViralLoad$log10VL <= Cens_VL) + 1*(Simulated_ViralLoad$obs %in% c(2,4))*(Simulated_ViralLoad$log10VL <= Cens_sgVL)
Simulated_ViralLoad$log10VL <- sapply(seq(1,nrow(Simulated_ViralLoad)),function(i)  max(Simulated_ViralLoad$log10VL[i],Cens_VL)*(Simulated_ViralLoad$obs[i] %in% c(1,3)) + max(Simulated_ViralLoad$log10VL[i],Cens_sgVL)*(Simulated_ViralLoad$obs[i] %in% c(2,4)))
# We only keep time points that are observed in real data
Kept_Simulated_ViralLoad <- Simulated_ViralLoad
Kept_Simulated_ViralLoad$log10VL <- sapply(seq(1,nrow(Kept_Simulated_ViralLoad)),function(i,data){ifelse((data$obs[i] %in% c(1,3))*(data$Time[i] %in% Time_VL) | (data$obs[i] %in% c(2,4))*(data$Time[i] %in% Time_sgVL),
                                                                                                         data$log10VL[i],NA)},data=Kept_Simulated_ViralLoad)
# Record of simulated viral load data and individual parameters
write.table(Kept_Simulated_ViralLoad,file=paste(Simulation_folder,"Simulated_data_VL_ECL_UnNoised.csv",sep="/"),row.names = FALSE,dec=".",sep=";")
write.table(Simulated_Parameters,file=paste(Simulation_folder,"Simulated_parameters_ViralLoad_beta_delta.csv",sep="/"),row.names=FALSE,dec=".",sep=";")
# ------------------------- #


# ------------------------- #
# PART III: Creation of datasets with VL and noised markers #### 
Simulated_noised_ECL <- read.table(file=paste(Simulation_folder,paste("Simulated_ECLRBD_ExternalData_Noised.csv",sep=""),sep="/"),header=TRUE,stringsAsFactors = FALSE,sep=";")
Simulated_noised_Marker2 <- read.table(file=paste(Simulation_folder,paste("Simulated_Marker2_ExternalData_Noised.csv",sep=""),sep="/"),header=TRUE,stringsAsFactors = FALSE,sep=";")
Simulated_noised_Marker3 <- read.table(file=paste(Simulation_folder,paste("Simulated_Marker3_ExternalData_Noised.csv",sep=""),sep="/"),header=TRUE,stringsAsFactors = FALSE,sep=";")
Kept_Simulated_ViralLoad <- read.table(file=paste(Simulation_folder,"Simulated_data_VL_ECL_UnNoised.csv",sep="/"),header = TRUE,stringsAsFactors = FALSE,dec=".",sep=";")

Combined_dataset <- Kept_Simulated_ViralLoad

# -- Addition of ECL in the dataset 
# - We calculate linear interpolation of the noised ECL data
LinearInterpolationParameters <- Linear_Interpolation_Parameters(Data=Simulated_noised_ECL,markers="ECL")
Simulated_noised_ECL$ECL.Slope <- sapply(seq(1,nrow(Simulated_noised_ECL)),function(i,param,data){
  tmp <- param$Slope[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=LinearInterpolationParameters,data=Simulated_noised_ECL)
Simulated_noised_ECL$ECL.Intercept <- sapply(seq(1,nrow(Simulated_noised_ECL)),function(i,param,data){
  tmp <- param$Intercept[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=LinearInterpolationParameters,data=Simulated_noised_ECL)
# We modify the columns related to ecl with the values of noised ecl
Combined_dataset$ECL <- NA ; Combined_dataset$ECL.Slope <- NA ; Combined_dataset$ECL.Intercept <- NA
# - We rebuild the linearly interpolated ECL noised dynamics at each timepoint of the observed viral load using parameters obtained only with ECL timepoints
Combined_dataset$ECL <- sapply(seq(1,nrow(Combined_dataset)),function(i,data,params){
  params_id <- params[which(params$AnimalID == data$AnimalID[i]),]
  if(data$Time[i] %in% params_id$Time){
    return(params_id$ECL[which(params_id$Time == data$Time[i])])
  }else{
    tmp_time <- max(params_id$Time[which(params_id$Time <= data$Time[i])],na.rm=TRUE)
    slope_time <- params_id$ECL.Slope[which(params_id$Time == tmp_time)]
    intercept_time <- params_id$ECL.Intercept[which(params_id$Time == tmp_time)]
    return(slope_time*data$Time[i]+intercept_time)
  }
},data=Combined_dataset,params=Simulated_noised_ECL)
# - We recalculate the parameters of the linear interpolation for all the timepoints of the viral load (needed afterwards for the model estimation in Monolix)
New_LinearInterpolationParameters <- Linear_Interpolation_Parameters(Data=Combined_dataset,markers="ECL")
Combined_dataset$ECL.Slope <- sapply(seq(1,nrow(Combined_dataset)),function(i,param,data){
  tmp <- param$Slope[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=New_LinearInterpolationParameters,data=Combined_dataset)
Combined_dataset$ECL.Intercept <- sapply(seq(1,nrow(Combined_dataset)),function(i,param,data){
  tmp <- param$Intercept[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=New_LinearInterpolationParameters,data=Combined_dataset)

# -- Addition of the 2nd marker in the dataset
LinearInterpolationParameters <- Linear_Interpolation_Parameters(Data=Simulated_noised_Marker2,markers="ECL")
Simulated_noised_Marker2$ECL.Slope <- sapply(seq(1,nrow(Simulated_noised_Marker2)),function(i,param,data){
  tmp <- param$Slope[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=LinearInterpolationParameters,data=Simulated_noised_Marker2)
Simulated_noised_Marker2$ECL.Intercept <- sapply(seq(1,nrow(Simulated_noised_Marker2)),function(i,param,data){
  tmp <- param$Intercept[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=LinearInterpolationParameters,data=Simulated_noised_Marker2)
# We modify the columns related to ecl with the values of noised ecl
Combined_dataset$Marker2 <- NA ; Combined_dataset$Marker2.Slope <- NA ; Combined_dataset$Marker2.Intercept <- NA
# - We rebuild the linearly interpolated ECL noised dynamics at each timepoint of the observed viral load using parameters obtained only with ECL timepoints
Combined_dataset$Marker2 <- sapply(seq(1,nrow(Combined_dataset)),function(i,data,params){
  params_id <- params[which(params$AnimalID == data$AnimalID[i]),]
  if(data$Time[i] %in% params_id$Time){
    return(params_id$ECL[which(params_id$Time == data$Time[i])])
  }else{
    tmp_time <- max(params_id$Time[which(params_id$Time <= data$Time[i])],na.rm=TRUE)
    slope_time <- params_id$ECL.Slope[which(params_id$Time == tmp_time)]
    intercept_time <- params_id$ECL.Intercept[which(params_id$Time == tmp_time)]
    return(slope_time*data$Time[i]+intercept_time)
  }
},data=Combined_dataset,params=Simulated_noised_Marker2)
# - We recalculate the parameters of the linear interpolation for all the timepoints of the viral load (needed afterwards for the model estimation in Monolix)
New_LinearInterpolationParameters <- Linear_Interpolation_Parameters(Data=Combined_dataset,markers="Marker2")
Combined_dataset$Marker2.Slope <- sapply(seq(1,nrow(Combined_dataset)),function(i,param,data){
  tmp <- param$Slope[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=New_LinearInterpolationParameters,data=Combined_dataset)
Combined_dataset$Marker2.Intercept <- sapply(seq(1,nrow(Combined_dataset)),function(i,param,data){
  tmp <- param$Intercept[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=New_LinearInterpolationParameters,data=Combined_dataset)

# -- Addition of the 2nd marker in the dataset
LinearInterpolationParameters <- Linear_Interpolation_Parameters(Data=Simulated_noised_Marker3,markers="ECL")
Simulated_noised_Marker3$ECL.Slope <- sapply(seq(1,nrow(Simulated_noised_Marker3)),function(i,param,data){
  tmp <- param$Slope[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=LinearInterpolationParameters,data=Simulated_noised_Marker3)
Simulated_noised_Marker3$ECL.Intercept <- sapply(seq(1,nrow(Simulated_noised_Marker3)),function(i,param,data){
  tmp <- param$Intercept[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=LinearInterpolationParameters,data=Simulated_noised_Marker3)
# We modify the columns related to ecl with the values of noised ecl
Combined_dataset$Marker3 <- NA ; Combined_dataset$Marker3.Slope <- NA ; Combined_dataset$Marker3.Intercept <- NA
# - We rebuild the linearly interpolated ECL noised dynamics at each timepoint of the observed viral load using parameters obtained only with ECL timepoints
Combined_dataset$Marker3 <- sapply(seq(1,nrow(Combined_dataset)),function(i,data,params){
  params_id <- params[which(params$AnimalID == data$AnimalID[i]),]
  if(data$Time[i] %in% params_id$Time){
    return(params_id$ECL[which(params_id$Time == data$Time[i])])
  }else{
    tmp_time <- max(params_id$Time[which(params_id$Time <= data$Time[i])],na.rm=TRUE)
    slope_time <- params_id$ECL.Slope[which(params_id$Time == tmp_time)]
    intercept_time <- params_id$ECL.Intercept[which(params_id$Time == tmp_time)]
    return(slope_time*data$Time[i]+intercept_time)
  }
},data=Combined_dataset,params=Simulated_noised_Marker3)
# - We recalculate the parameters of the linear interpolation for all the timepoints of the viral load (needed afterwards for the model estimation in Monolix)
New_LinearInterpolationParameters <- Linear_Interpolation_Parameters(Data=Combined_dataset,markers="Marker3")
Combined_dataset$Marker3.Slope <- sapply(seq(1,nrow(Combined_dataset)),function(i,param,data){
  tmp <- param$Slope[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=New_LinearInterpolationParameters,data=Combined_dataset)
Combined_dataset$Marker3.Intercept <- sapply(seq(1,nrow(Combined_dataset)),function(i,param,data){
  tmp <- param$Intercept[which(param$Time == data$Time[i] & param$AnimalID == data$AnimalID[i])]
  return(tmp)
},param=New_LinearInterpolationParameters,data=Combined_dataset)

# Rearragement of the dataset
Combined_dataset <- Combined_dataset[,c("AnimalID","Time","log10VL","obs","Group","type","organ","Weight","censored",paste("ECL",c("",".Slope",".Intercept"),sep=""),paste("Marker2",c("",".Slope",".Intercept"),sep=""),paste("Marker3",c("",".Slope",".Intercept"),sep=""))]

# We save the FINAL dataset
write.table(Combined_dataset,file=paste(Simulation_folder,paste("Simulated_data_ViralLoad_ImmunoMarkers.csv",sep=""),sep="/"),row.names = FALSE,dec=".",sep=";")
# ------------------------- #