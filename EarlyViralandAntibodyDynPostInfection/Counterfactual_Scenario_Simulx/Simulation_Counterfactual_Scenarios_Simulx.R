# --------------- # 
# DESCRIPTION: Simulation of counterfactual scenarios using the estimated Joint model and Simulx
#
# Author: Marie Alexandre 
# R version: 4.2.1
# Simulx Version: 2023R1


# Open the document outline (Ctrl + shift + O) to see the structure of the document
# --------------- # 



rm(list=ls())
'%notin%' <- Negate('%in%')





# --- FUNCTIONS ----
# -- General functions --
Plot_colors_Groups <- function(){
  colors <-  c(rgb(128,128,128,maxColorValue = 255),rgb(181,219,69,maxColorValue = 255),rgb(69,181,219,maxColorValue = 255),
               rgb(219,69,181,maxColorValue = 255),rgb(219,107,69,maxColorValue = 255),rgb(80,47,168,maxColorValue = 255))
  names(colors) <- c("Naive","CD40.RBDv - Naive","Conv","CD40.RBDv - Conv","mRNA - Conv","CD40.PanCoV - Conv")
  return(colors)
}
exp_bold <- function(lab) {
  new_lab <- log10(lab)
  do.call(
    expression,
    lapply(paste(new_lab), function(x) bquote(bold("10"^.(x))))
  )
}

Distribution_Function <- function(.data,.variables){
  Distribution <- ddply(.data=.data,.variables = .variables,summarize,
                        Mean=mean(Value,na.rm=TRUE),
                        Median=median(Value,na.rm=TRUE),
                        ICMIN=quantile(Value,na.rm=TRUE,probs=c(0.025)),
                        Q1=quantile(Value,na.rm=TRUE,probs=c(0.25)),
                        Q3=quantile(Value,na.rm=TRUE,probs=c(0.75)),
                        ICMAX=quantile(Value,na.rm=TRUE,probs=c(0.975)))
  return(Distribution)
}

# -- Function descriptors --
Peak_VL_dynamics <- function(time,value){
  data <- data.frame(time=time,value=value,row.names = NULL)
  data <- data[order(data$time),]
  
  return(max(data$value))
}
TimePeak_VL_dynamics <- function(time,value){
  data <- data.frame(time=time,value=value,row.names = NULL)
  data <- data[order(data$time),]
  
  return(data$time[which.max(data$value)])
}
AUC_VL_Dynamics <- function(time,value){
  data <- data.frame(time=time,value=value,row.names = NULL)
  data <- data[order(data$time),]
  auc_value <- AUC(x=data$time,y=data$value,method="trapezoid")
  
  return(auc_value)
}
Duration_Clearance_Stage <- function(time,value,lod){
  # time interval between the peak and the first undetectable vl
  data <- data.frame(time=time,value=value,row.names = NULL)
  data <- data[order(data$time),]
  data$censored <- 1*(data$value <= lod)
  
  time_peak <- TimePeak_VL_dynamics(time=time,value=value)
  duration_clearance <- min(subset(data,time>=time_peak & censored == 1)$time) - time_peak
  
  return(duration_clearance)
}
Duration_Acute_Stage <- function(time,value,lod){
  # time interval between the first and the last detectable vl
  data <- data.frame(time=time,value=value,row.names = NULL)
  data <- data[order(data$time),]
  data$censored <- 1*(data$value <= lod)
  
  if(sum(data$censored) == nrow(data)){
    duration_acute <- 0
  }else if(data$censored[nrow(data)] == 0){
    duration_acute <- Inf
  }else{
    duration_acute <- diff(range(subset(data,censored==0)$time))
  }
  
  return(duration_acute)
}

Plot_Real_Vs_Counterfactual_Dynamics <- function(data,facet_labels){
  
  Colors <- Plot_colors_Groups()
  data$Group <- factor(data$Group,levels = names(Colors))
  strip_settings <- strip_themed(background_x = elem_list_rect(fill=alpha(Colors,alpha = 0.25),color=alpha(Colors,alpha = 1)),
                                 by_layer_x = FALSE,
                                 background_y = element_rect(fill="white",color="white"),
                                 text_y = element_text(color="black",size=11,face="bold"))
  
  
  plot <-  ggplot(data=data) + 
    # Confidence Interval of counterfactual dynamics
    geom_ribbon(data=subset(data,Type == "Counterfactual"),
                aes(x=time,ymin=ICMIN,ymax=ICMAX,fill=as.factor(Group),color=as.factor(Group)),linetype="solid",linewidth=0.80,alpha=0.3) +
    geom_ribbon(data=subset(data,Type == "Counterfactual"),
                aes(x=time,ymin=ICMIN,ymax=ICMAX,color=as.factor(Group)),linetype="solid",fill=NA,linewidth=0.80,alpha=1) +
    
    # Confidence Interval of real dynamics
    geom_ribbon(data=subset(data,Type == "Real"),
                aes(x=time,ymin=ICMIN,ymax=ICMAX),color="black",linetype="dashed",fill=NA,linewidth=0.80,alpha=0.3) +
    geom_ribbon(data=subset(data,Type == "Real"),
                aes(x=time,ymin=ICMIN,ymax=ICMAX),color="black",linetype="dashed",fill=NA,linewidth=0.80,alpha=1) +
    
    # Mean of counterfactual dynamics
    geom_line(data=subset(data,Type == "Counterfactual"),
              aes(x=time,y=Median,color=as.factor(Group),linetype=as.factor(Type)),linewidth=1.2) +
    # Mean of real dynamics
    geom_line(data=subset(data,Type == "Real"),
              aes(x=time,y=Median,linetype=as.factor(Type)),color="black",linewidth=1.2) +
    
    facet_nested(ObservationType~Group,scales = "free_y",strip=strip_settings,switch="y",labeller=labeller(ObservationType=facet_labels)) +
    
    scale_color_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE) +
    scale_fill_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE)  +
    scale_linetype_manual(name="Type of dynamics",breaks = c("Counterfactual","Real"),values = c("solid","dashed"),labels = c("Counterfactual dynamics","Real dynamics")) + 
    guides(fill = FALSE, color=FALSE,linetype=guide_legend(nrow=1)) +
    
    theme_bw() +
    theme(axis.line = element_line(color="black",size=1.0),
          axis.text =  element_text(color="black",size=10,face="bold"),
          axis.ticks = element_line(color="black",size=1),
          axis.ticks.length = unit(0.1,units = "cm"),
          axis.title = element_text(color="black",size=11,face="bold"),
          legend.title = element_text(face = "bold"),
          legend.text = element_text(color="black",size=11),
          legend.key = element_rect(color="white",fill="white",size=0.5),
          legend.key.width = unit(units = "cm",x = 1.5),
          legend.key.height = unit(units = "cm",x = 0.1),
          legend.position = "top",
          legend.direction = "horizontal",
          panel.spacing.x = unit(0,units="cm"),
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size = 10,face="bold"))  +
    theme(strip.placement.y = "outside",
          strip.text.y = element_text(margin = margin(unit = "cm",r = 0,l=0)),
          axis.title.y = element_blank())
  
  
  return(plot)
}
# ------------------ #



# --- LIBRAIRIES ----
library(lixoftConnectors,lib.loc = "C:/Users/marie/AppData/Local/R/win-library/4.4/Monolix2023")
library(RsSimulx)
# 
library(DescTools)  # Library for the function AUC 

# # Libraries for visual representations
library(ggplot2)
library(lemon)
library(scales)   # Library for scientific notation in y-axis
library(ggh4x)

library(reshape2) # Library to reshape dataset
library(plyr)  # Library used to calculate data distribution
# ------------------ #







# -- Dataset --
Project_Folder <- "EarlyViralandAntibodyDynPostInfection"

Data_Folder <- paste(Project_Folder,"Simulated_data",sep="/")
Data_file <- "Simulated_dataset_VL_BAb_NAb_Delta_PostExpo.txt"
# Load of dataset
ObservedData <- read.table(file = paste(Data_Folder,Data_file,sep="/"),header=TRUE,sep="\t",dec=".")
ObservedData$MeasureType[which(ObservedData$obsid == 5)] <- "BAb"
ObservedData$MeasureType[which(ObservedData$obsid == 6)] <- "NAb"
NHPs <- unique(ObservedData[,c("Group","SubjectID","Inoc","Weight")])

LOD_gRNA <- 476 # in copies/mL
LOD_sgRNA <- 749 # in copies/mL

# -- Monolix Model -- 
MonolixModel_Folder <- paste(Project_Folder,"Model_Estimation_Monolix",sep="/")
MonolixModel_name <- "Final_JointModel_NHP_DeltaInfection_MONOLIX"
Monolix_Project <- paste(MonolixModel_Folder,paste(MonolixModel_name,"mlxtran",sep="."),sep="/")

Results_Folder <- paste(Project_Folder,"Counterfactual_Scenario_Simulx",sep="/")

# -- Figures ---
Figure_Folder <- "Figures"
# Uncomment to save figures 
# dir.create(paste(Results_Folder,Figure_Folder,sep="/"),recursive = TRUE)

# -- Simulation Code --
MlxtranModel_CounterfactualScenario_Ablevel_file <- paste(Results_Folder,"MlxtranCode_Model_gRNA_sgRNA_BAb_NAb_PostInfection_CounterFactualScenario_Ablevel_SIMULATION.txt",sep="/")
MlxtranModel_CounterfactualScenario_Nodepl_file <- paste(Results_Folder,"MlxtranCode_Model_gRNA_sgRNA_BAb_NAb_PostInfection_CounterFactualScenario_Nodepletion_SIMULATION.txt",sep="/")
MlxtranModel_CounterfactualScenario_NoNeut_file <- paste(Results_Folder,"MlxtranCode_Model_gRNA_sgRNA_BAb_NAb_PostInfection_CounterFactualScenario_NoNeutralization_SIMULATION.txt",sep="/")



# -- Parameters of simulation ---
Nb_pop <- 500     # Number of population parameters for median dynamics
Seed <- 123456    # Seed used to simulated parameters for reproductibility of the work
Time_simulation <- seq(0,30,by=0.25)
Model_Outputs <- list(list(name="Obs_NAb",time=Time_simulation),
                      list(name="gRNA_N",time=Time_simulation),
                      list(name="sgRNA_N",time=Time_simulation),
                      list(name="TN_T0_Percent",time=Time_simulation))

NHPs_Profile <- data.frame(Group = unique(NHPs$Group),Inoc=median(NHPs$Inoc),Weight=median(NHPs$Weight))
NHPs_Profile <- cbind(id = seq(1,nrow(NHPs_Profile)),NHPs_Profile)

Covariates <- NHPs_Profile[,c("id","Group")]
Regressors <- merge(NHPs_Profile,expand.grid(id=unique(NHPs_Profile$id),time=Time_simulation))
Regressors <- Regressors[order(Regressors$time,Regressors$id),]
# ------------ #






# ............................----
# -- SIMULATION OF POPULATION PARAMETERS  ----
# ............................----
initializeLixoftConnectors(software = "simulx",path = "C:/ProgramData/Lixoft/MonolixSuite2023R1",force = T)
Simulated_Population_Parameters <- simpopmlx(n=Nb_pop,project = Monolix_Project,seed = Seed)
# Removal of inter-individual variability and observation uncertainty
Simulated_Population_Parameters[,which(names(Simulated_Population_Parameters) %in% c(paste("omega",c("delta_N","S0","theta_BAb","eta"),sep="_"),
                                                                                     paste("a",seq(1,6),sep="")))] <- 0
# ------------ #



# ............................----
# -- SIMULATION OF REAL DYNAMICS ----
# ............................----
MlxtranModel_Realmodel_file <- paste(MonolixModel_Folder,"MlxtranCode_NHP_JointModel_VL_Ab_PostInfection_SIMULX.txt",sep="/")

### > Simulation of dynamics ----
Simulated_Real_Dynamics <- simulx(model = MlxtranModel_Realmodel_file,
                                  parameter = cbind(Simulated_Population_Parameters),
                                  output = Model_Outputs,regressor = Regressors,
                                  covariate = Covariates,
                                  settings = list(sharedIds="population"))[sapply(Model_Outputs, function(l) l$name)]

Simulated_Real_Dynamics <- Reduce(function(x,y) merge(x,y,all=TRUE), Simulated_Real_Dynamics)
Simulated_Real_Dynamics <- merge(Covariates,Simulated_Real_Dynamics)[-1]
names(Simulated_Real_Dynamics)[which(names(Simulated_Real_Dynamics) == "rep")] <- "pop"
Simulated_Real_Dynamics$Obs_NAb <- 10^Simulated_Real_Dynamics$Obs_NAb
Simulated_Real_Dynamics$gRNA_N <- 10^Simulated_Real_Dynamics$gRNA_N
Simulated_Real_Dynamics$sgRNA_N <- 10^Simulated_Real_Dynamics$sgRNA_N

### > Estimation of VL descriptors ----
Descriptors_Real_gRNA_Dynamics <- ddply(.data=Simulated_Real_Dynamics,.variables = .(pop,Group),summarize,
                                        Peak=Peak_VL_dynamics(time=time,value=gRNA_N),
                                        TimePeak=TimePeak_VL_dynamics(time=time,value=gRNA_N),
                                        AUC=AUC_VL_Dynamics(time=time,value=gRNA_N),
                                        Clearance=Duration_Clearance_Stage(time=time,value=gRNA_N,lod=LOD_gRNA),
                                        Acute=Duration_Acute_Stage(time=time,value=gRNA_N,lod=LOD_gRNA))
Descriptors_Real_gRNA_Dynamics <- cbind(VLType="gRNA",Descriptors_Real_gRNA_Dynamics)

Descriptors_Real_sgRNA_Dynamics <- ddply(.data=Simulated_Real_Dynamics,.variables = .(pop,Group),summarize,
                                         Peak=Peak_VL_dynamics(time=time,value=sgRNA_N),
                                         TimePeak=TimePeak_VL_dynamics(time=time,value=sgRNA_N),
                                         AUC=AUC_VL_Dynamics(time=time,value=sgRNA_N),
                                         Clearance=Duration_Clearance_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA),
                                         Acute=Duration_Acute_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA))
Descriptors_Real_sgRNA_Dynamics <- cbind(VLType="sgRNA",Descriptors_Real_sgRNA_Dynamics)
# ------------ #




# ............................----
# -- COUNTERFACTUAL SCENARIO AB LEVEL AT EXPOSURE ----
# ............................----
# Scenario: What if NAb at exposure was equal to the same value in all groups ?


## Scenario LownAb: Low level of antibodies ----
Counterfactual_NAb0_1a <- 0.90  # Protective threshold estimated in PanCoV group


### > Simulation of dynamics ---- 
Simulated_Scenario1a_Dynamics <- simulx(model = MlxtranModel_CounterfactualScenario_Ablevel_file,
                                        parameter = cbind(Simulated_Population_Parameters,NAb_thres = Counterfactual_NAb0_1a),
                                        output = Model_Outputs,regressor = Regressors,
                                        covariate = Covariates,
                                        settings = list(sharedIds="population"))[sapply(Model_Outputs, function(l) l$name)]

Simulated_Scenario1a_Dynamics <- Reduce(function(x,y) merge(x,y,all=TRUE), Simulated_Scenario1a_Dynamics)
Simulated_Scenario1a_Dynamics <- merge(Covariates,Simulated_Scenario1a_Dynamics)[-1]
names(Simulated_Scenario1a_Dynamics)[which(names(Simulated_Scenario1a_Dynamics) == "rep")] <- "pop"
Simulated_Scenario1a_Dynamics$Obs_NAb <- 10^Simulated_Scenario1a_Dynamics$Obs_NAb
Simulated_Scenario1a_Dynamics$gRNA_N <- 10^Simulated_Scenario1a_Dynamics$gRNA_N
Simulated_Scenario1a_Dynamics$sgRNA_N <- 10^Simulated_Scenario1a_Dynamics$sgRNA_N


### > Estimation of VL descriptors ----
Descriptors_Scenario1a_gRNA_Dynamics <- ddply(.data=Simulated_Scenario1a_Dynamics,.variables = .(pop,Group),summarize,
                                              Peak=Peak_VL_dynamics(time=time,value=gRNA_N),
                                              TimePeak=TimePeak_VL_dynamics(time=time,value=gRNA_N),
                                              AUC=AUC_VL_Dynamics(time=time,value=gRNA_N),
                                              Clearance=Duration_Clearance_Stage(time=time,value=gRNA_N,lod=LOD_gRNA),
                                              Acute=Duration_Acute_Stage(time=time,value=gRNA_N,lod=LOD_gRNA))
Descriptors_Scenario1a_gRNA_Dynamics <- cbind(VLType="gRNA",Descriptors_Scenario1a_gRNA_Dynamics)

Descriptors_Scenario1a_sgRNA_Dynamics <- ddply(.data=Simulated_Scenario1a_Dynamics,.variables = .(pop,Group),summarize,
                                               Peak=Peak_VL_dynamics(time=time,value=sgRNA_N),
                                               TimePeak=TimePeak_VL_dynamics(time=time,value=sgRNA_N),
                                               AUC=AUC_VL_Dynamics(time=time,value=sgRNA_N),
                                               Clearance=Duration_Clearance_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA),
                                               Acute=Duration_Acute_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA))
Descriptors_Scenario1a_sgRNA_Dynamics <- cbind(VLType="sgRNA",Descriptors_Scenario1a_sgRNA_Dynamics)



### > Plot of dynamics ----
# Reshape of the dataset
Reshaped_Simulated_Scenario1a_Dynamics <- melt(data=Simulated_Scenario1a_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")
Reshaped_Simulated_Real_Dynamics <- melt(data=Simulated_Real_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")
Merged_Simulated_Dynamics_Scenario1a <- rbind(cbind(Type="Counterfactual",Reshaped_Simulated_Scenario1a_Dynamics),
                                              cbind(Type="Real",Reshaped_Simulated_Real_Dynamics))


# Calculation of distributions
Distribution_Dynamics_Scenario1a <- Distribution_Function(.data=Merged_Simulated_Dynamics_Scenario1a,.variables = .(Type,Group,ObservationType,time))

facet_labels <- c("ACE2/RBD bind. inhibit° \n (AU/mL)","Genomic RNA \n (copies/mL)","Subgenomic RNA \n (copies/mL)","Target cells \n  (T/T0, %)")
names(facet_labels) <- c("Obs_NAb","gRNA_N","sgRNA_N","TN_T0_Percent")

LOD_values <- data.frame(ObservationType=names(facet_labels),Value=c(0.1,LOD_gRNA,LOD_sgRNA,NA))
LOD_values$ObservationType <- factor(LOD_values$ObservationType,levels=names(facet_labels))

# --- Viral load dynamics (Comparison with real dynamics)
Global_Plot_Scenario1a <- Plot_Real_Vs_Counterfactual_Dynamics(data=subset(Distribution_Dynamics_Scenario1a,ObservationType != "Obs_NAb"),facet_labels = facet_labels) + 
  geom_hline(data=subset(LOD_values,ObservationType != "Obs_NAb"),
             aes(yintercept = Value),color="darkred",linetype="solid",linewidth=0.8,alpha=0.5) +
  scale_x_continuous(name="Time post-exposure (days)",breaks = seq(0,30,by=5))  +
  facetted_pos_scales(y = list(
    ObservationType == "Obs_NAb" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(2E-2,0.5E3)),
    ObservationType == "gRNA_N" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(1E-2,1E9)),
    ObservationType == "sgRNA_N" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(1E-3,5E7)),
    ObservationType == "TN_T0_Percent" ~ scale_y_continuous(breaks = seq(0,100,by=20), limits=c(-1,101)))) + 
  coord_cartesian(xlim = c(-0.5,29.5),expand=0)
# ggsave(Global_Plot_Scenario1a,filename=paste(Results_Folder,Figure_Folder,"ViralPlots_Counterfactual_Scenario_LownAb.png",sep="/"),height=7,width=10,dpi=300)


# --- Neutralizing antibody dynamics
Colors <- Plot_colors_Groups()
NAb_Plot_Scenario1a <- ggplot(subset(Distribution_Dynamics_Scenario1a,Type == "Counterfactual" & ObservationType == "Obs_NAb")) + 
  geom_ribbon(aes(x=time,ymin=ICMIN,ymax=ICMAX,fill=as.factor(Group),color=as.factor(Group)),linetype="dotted",linewidth=0.30,alpha=0.3) +
  geom_ribbon(aes(x=time,ymin=ICMIN,ymax=ICMAX,color=as.factor(Group)),linetype="dotted",fill=NA,linewidth=0.30,alpha=1) +
  geom_line(aes(x=time,y=Mean,color=as.factor(Group)),linewidth=1.0) + 
  
  scale_color_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE) +
  scale_fill_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE)  +
  guides(color=guide_legend(nrow=6)) + 
  
  scale_x_continuous(name="Time post-exposure (days)",breaks = seq(0,30,by=5))  +
  scale_y_log10(breaks = c(0.5,1,2,3,5,10,20,30,50,100)) + 
  ylab("ACE2/RBD binding inhibit° (AU/mL)") + 
  coord_cartesian(xlim = c(-0.5,25.5),ylim=c(0.5,0.5E2),expand=0)  +
  ggtitle(expression(bold('NAb'[0]*'= 0.90'*' AU/mL'))) + 
  
  theme_bw() +
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=11,face="bold"),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(color="black",face="bold",size=10),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        legend.key.width = unit(units = "cm",x = 1.25),
        legend.key.height = unit(units = "cm",x = 0.1),
        legend.position = "inside",
        legend.position.inside = c(0.78,0.3),
        legend.title.position = "top",
        legend.direction = "horizontal",
        panel.spacing.x = unit(0,units="cm"),
        panel.grid.major = element_line(colour = "gray90"))
# ggsave(NAb_Plot_Scenario1a,filename=paste(Results_Folder,Figure_Folder,"NAbPlot_Counterfactual_Scenario_LownAb.png",sep="/"),height=3.5,width=6,dpi=300)




## Scenario MidnAb: Medium level of antibodies ----
Counterfactual_NAb0_1b <- median(subset(ObservedData,Group=="CD40.PanCoV - Conv")$NAb0)

### > Simulation of dynamics ---- 
Simulated_Scenario1b_Dynamics <- simulx(model = MlxtranModel_CounterfactualScenario_Ablevel_file,
                                        parameter = cbind(Simulated_Population_Parameters,NAb_thres = Counterfactual_NAb0_1b),
                                        output = Model_Outputs,regressor = Regressors,
                                        covariate = Covariates,
                                        settings = list(sharedIds="population"))[sapply(Model_Outputs, function(l) l$name)]

Simulated_Scenario1b_Dynamics <- Reduce(function(x,y) merge(x,y,all=TRUE), Simulated_Scenario1b_Dynamics)
Simulated_Scenario1b_Dynamics <- merge(Covariates,Simulated_Scenario1b_Dynamics)[-1]
names(Simulated_Scenario1b_Dynamics)[which(names(Simulated_Scenario1b_Dynamics) == "rep")] <- "pop"
Simulated_Scenario1b_Dynamics$Obs_NAb <- 10^Simulated_Scenario1b_Dynamics$Obs_NAb
Simulated_Scenario1b_Dynamics$gRNA_N <- 10^Simulated_Scenario1b_Dynamics$gRNA_N
Simulated_Scenario1b_Dynamics$sgRNA_N <- 10^Simulated_Scenario1b_Dynamics$sgRNA_N


### > Estimation of VL descriptors ----
Descriptors_Scenario1b_gRNA_Dynamics <- ddply(.data=Simulated_Scenario1b_Dynamics,.variables = .(pop,Group),summarize,
                                              Peak=Peak_VL_dynamics(time=time,value=gRNA_N),
                                              TimePeak=TimePeak_VL_dynamics(time=time,value=gRNA_N),
                                              AUC=AUC_VL_Dynamics(time=time,value=gRNA_N),
                                              Clearance=Duration_Clearance_Stage(time=time,value=gRNA_N,lod=LOD_gRNA),
                                              Acute=Duration_Acute_Stage(time=time,value=gRNA_N,lod=LOD_gRNA))
Descriptors_Scenario1b_gRNA_Dynamics <- cbind(VLType="gRNA",Descriptors_Scenario1b_gRNA_Dynamics)

Descriptors_Scenario1b_sgRNA_Dynamics <- ddply(.data=Simulated_Scenario1b_Dynamics,.variables = .(pop,Group),summarize,
                                               Peak=Peak_VL_dynamics(time=time,value=sgRNA_N),
                                               TimePeak=TimePeak_VL_dynamics(time=time,value=sgRNA_N),
                                               AUC=AUC_VL_Dynamics(time=time,value=sgRNA_N),
                                               Clearance=Duration_Clearance_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA),
                                               Acute=Duration_Acute_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA))
Descriptors_Scenario1b_sgRNA_Dynamics <- cbind(VLType="sgRNA",Descriptors_Scenario1b_sgRNA_Dynamics)



### > Plot of dynamics ----
# Reshape of the dataset
Reshaped_Simulated_Scenario1b_Dynamics <- melt(data=Simulated_Scenario1b_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")
Reshaped_Simulated_Real_Dynamics <- melt(data=Simulated_Real_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")

Merged_Simulated_Dynamics_Scenario1b <- rbind(cbind(Type="Counterfactual",Reshaped_Simulated_Scenario1b_Dynamics),
                                              cbind(Type="Real",Reshaped_Simulated_Real_Dynamics))
# Calculation of distributions
Distribution_Dynamics_Scenario1b <- Distribution_Function(.data=Merged_Simulated_Dynamics_Scenario1b,.variables = .(Type,Group,ObservationType,time))
Distribution_Dynamics_Scenario1b$ObservationType <- factor(Distribution_Dynamics_Scenario1b$ObservationType,levels=c("Obs_NAb","gRNA_N","sgRNA_N","TN_T0_Percent"))


facet_labels <- c("ACE2/RBD bind. inhibit° \n (AU/mL)","Genomic RNA \n (copies/mL)","Subgenomic RNA \n (copies/mL)","Target cells \n  (T/T0, %)")
names(facet_labels) <- c("Obs_NAb","gRNA_N","sgRNA_N","TN_T0_Percent")

LOD_values <- data.frame(ObservationType=names(facet_labels),Value=c(0.1,LOD_gRNA,LOD_sgRNA,NA))
LOD_values$ObservationType <- factor(LOD_values$ObservationType,levels=names(facet_labels))

Global_Plot_Scenario1b <- Plot_Real_Vs_Counterfactual_Dynamics(data=subset(Distribution_Dynamics_Scenario1b,ObservationType != "Obs_NAb"),facet_labels = facet_labels) + 
  geom_hline(data=subset(LOD_values,ObservationType != "Obs_NAb"),
             aes(yintercept = Value),color="darkred",linetype="solid",linewidth=0.8,alpha=0.5) +
  scale_x_continuous(name="Time post-exposure (days)",breaks = seq(0,30,by=5))  +
  facetted_pos_scales(y = list(
    ObservationType == "Obs_NAb" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(2E-2,0.5E3)),
    ObservationType == "gRNA_N" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(1E-2,1E9)),
    ObservationType == "sgRNA_N" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(1E-3,5E7)),
    ObservationType == "TN_T0_Percent" ~ scale_y_continuous(breaks = seq(0,100,by=20), limits=c(-1,101)))) + 
  coord_cartesian(xlim = c(-0.5,29.5),expand=0) 
# ggsave(Global_Plot_Scenario1b,filename=paste(Results_Folder,Figure_Folder,"ViralPlots_Counterfactual_Scenario_MidnAb.png",sep="/"),height=7,width=10,dpi=300)




## Scenario HinAb: High level of antibodies ----
Counterfactual_NAb0_1c <- max(subset(ObservedData,Group=="CD40.PanCoV - Conv")$NAb0)

### > Simulation of dynamics ---- 
Simulated_Scenario1c_Dynamics <- simulx(model = MlxtranModel_CounterfactualScenario_Ablevel_file,
                                        parameter = cbind(Simulated_Population_Parameters,NAb_thres = Counterfactual_NAb0_1c),
                                        output = Model_Outputs,regressor = Regressors,
                                        covariate = Covariates,
                                        settings = list(sharedIds="population"))[sapply(Model_Outputs, function(l) l$name)]

Simulated_Scenario1c_Dynamics <- Reduce(function(x,y) merge(x,y,all=TRUE), Simulated_Scenario1c_Dynamics)
Simulated_Scenario1c_Dynamics <- merge(Covariates,Simulated_Scenario1c_Dynamics)[-1]
names(Simulated_Scenario1c_Dynamics)[which(names(Simulated_Scenario1c_Dynamics) == "rep")] <- "pop"
Simulated_Scenario1c_Dynamics$Obs_NAb <- 10^Simulated_Scenario1c_Dynamics$Obs_NAb
Simulated_Scenario1c_Dynamics$gRNA_N <- 10^Simulated_Scenario1c_Dynamics$gRNA_N
Simulated_Scenario1c_Dynamics$sgRNA_N <- 10^Simulated_Scenario1c_Dynamics$sgRNA_N


### > Estimation of VL descriptors ----
Descriptors_Scenario1c_gRNA_Dynamics <- ddply(.data=Simulated_Scenario1c_Dynamics,.variables = .(pop,Group),summarize,
                                              Peak=Peak_VL_dynamics(time=time,value=gRNA_N),
                                              TimePeak=TimePeak_VL_dynamics(time=time,value=gRNA_N),
                                              AUC=AUC_VL_Dynamics(time=time,value=gRNA_N),
                                              Clearance=Duration_Clearance_Stage(time=time,value=gRNA_N,lod=LOD_gRNA),
                                              Acute=Duration_Acute_Stage(time=time,value=gRNA_N,lod=LOD_gRNA))
Descriptors_Scenario1c_gRNA_Dynamics <- cbind(VLType="gRNA",Descriptors_Scenario1c_gRNA_Dynamics)

Descriptors_Scenario1c_sgRNA_Dynamics <- ddply(.data=Simulated_Scenario1c_Dynamics,.variables = .(pop,Group),summarize,
                                               Peak=Peak_VL_dynamics(time=time,value=sgRNA_N),
                                               TimePeak=TimePeak_VL_dynamics(time=time,value=sgRNA_N),
                                               AUC=AUC_VL_Dynamics(time=time,value=sgRNA_N),
                                               Clearance=Duration_Clearance_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA),
                                               Acute=Duration_Acute_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA))
Descriptors_Scenario1c_sgRNA_Dynamics <- cbind(VLType="sgRNA",Descriptors_Scenario1c_sgRNA_Dynamics)


### > Plot of dynamics ----
# Reshape of the dataset
Reshaped_Simulated_Scenario1c_Dynamics <- melt(data=Simulated_Scenario1c_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")
Reshaped_Simulated_Real_Dynamics <- melt(data=Simulated_Real_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")

Merged_Simulated_Dynamics_Scenario1c <- rbind(cbind(Type="Counterfactual",Reshaped_Simulated_Scenario1c_Dynamics),
                                              cbind(Type="Real",Reshaped_Simulated_Real_Dynamics))
# Calculation of distributions
Distribution_Dynamics_Scenario1c <- Distribution_Function(.data=Merged_Simulated_Dynamics_Scenario1c,.variables = .(Type,Group,ObservationType,time))
Distribution_Dynamics_Scenario1c$ObservationType <- factor(Distribution_Dynamics_Scenario1c$ObservationType,levels=c("Obs_NAb","gRNA_N","sgRNA_N","TN_T0_Percent"))


facet_labels <- c("ACE2/RBD bind. inhibit° \n (AU/mL)","Genomic RNA \n (copies/mL)","Subgenomic RNA \n (copies/mL)","Target cells \n  (T/T0, %)")
names(facet_labels) <- c("Obs_NAb","gRNA_N","sgRNA_N","TN_T0_Percent")

LOD_values <- data.frame(ObservationType=names(facet_labels),Value=c(0.1,LOD_gRNA,LOD_sgRNA,NA))
LOD_values$ObservationType <- factor(LOD_values$ObservationType,levels=names(facet_labels))

Global_Plot_Scenario1c <- Plot_Real_Vs_Counterfactual_Dynamics(data=subset(Distribution_Dynamics_Scenario1c,ObservationType != "Obs_NAb"),facet_labels = facet_labels) + 
  geom_hline(data=subset(LOD_values,ObservationType != "Obs_NAb"),
             aes(yintercept = Value),color="darkred",linetype="solid",linewidth=0.8,alpha=0.5) +
  scale_x_continuous(name="Time post-exposure (days)",breaks = seq(0,30,by=5))  +
  facetted_pos_scales(y = list(
    ObservationType == "Obs_NAb" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(2E-2,0.5E3)),
    ObservationType == "gRNA_N" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(1E-2,1E9)),
    ObservationType == "sgRNA_N" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(1E-3,5E7)),
    ObservationType == "TN_T0_Percent" ~ scale_y_continuous(breaks = seq(0,100,by=20), limits=c(-1,101)))) + 
  coord_cartesian(xlim = c(-0.5,29.5),expand=0) 
# ggsave(Global_Plot_Scenario1c,filename=paste(Results_Folder,Figure_Folder,"ViralPlots_Counterfactual_Scenario_HinAb.png",sep="/"),height=7,width=10,dpi=300)
# ------------ #




# ............................----
# -- COUNTERFACTUAL SCENARIO NO TARGET-CELL DEPLETION ----
# ............................----
# Scenario: What if no target-cell depletion (T = T0 for all t) ?

### > Simulation of dynamics ---- 
Simulated_Scenario2_Dynamics <- simulx(model = MlxtranModel_CounterfactualScenario_Nodepl_file,
                                       parameter = cbind(Simulated_Population_Parameters),
                                       output = Model_Outputs,regressor = Regressors,
                                       covariate = Covariates,
                                       settings = list(sharedIds="population"))[sapply(Model_Outputs, function(l) l$name)]

Simulated_Scenario2_Dynamics <- Reduce(function(x,y) merge(x,y,all=TRUE), Simulated_Scenario2_Dynamics)
Simulated_Scenario2_Dynamics <- merge(Covariates,Simulated_Scenario2_Dynamics)[-1]
names(Simulated_Scenario2_Dynamics)[which(names(Simulated_Scenario2_Dynamics) == "rep")] <- "pop"
Simulated_Scenario2_Dynamics$Obs_NAb <- 10^Simulated_Scenario2_Dynamics$Obs_NAb
Simulated_Scenario2_Dynamics$gRNA_N <- 10^Simulated_Scenario2_Dynamics$gRNA_N
Simulated_Scenario2_Dynamics$sgRNA_N <- 10^Simulated_Scenario2_Dynamics$sgRNA_N


### > Estimation of VL descriptors ----
Descriptors_Scenario2_gRNA_Dynamics <- ddply(.data=Simulated_Scenario2_Dynamics,.variables = .(pop,Group),summarize,
                                             Peak=Peak_VL_dynamics(time=time,value=gRNA_N),
                                             TimePeak=TimePeak_VL_dynamics(time=time,value=gRNA_N),
                                             AUC=AUC_VL_Dynamics(time=time,value=gRNA_N),
                                             Clearance=Duration_Clearance_Stage(time=time,value=gRNA_N,lod=LOD_gRNA),
                                             Acute=Duration_Acute_Stage(time=time,value=gRNA_N,lod=LOD_gRNA))
Descriptors_Scenario2_gRNA_Dynamics <- cbind(VLType="gRNA",Descriptors_Scenario2_gRNA_Dynamics)

Descriptors_Scenario2_sgRNA_Dynamics <- ddply(.data=Simulated_Scenario2_Dynamics,.variables = .(pop,Group),summarize,
                                              Peak=Peak_VL_dynamics(time=time,value=sgRNA_N),
                                              TimePeak=TimePeak_VL_dynamics(time=time,value=sgRNA_N),
                                              AUC=AUC_VL_Dynamics(time=time,value=sgRNA_N),
                                              Clearance=Duration_Clearance_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA),
                                              Acute=Duration_Acute_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA))
Descriptors_Scenario2_sgRNA_Dynamics <- cbind(VLType="sgRNA",Descriptors_Scenario2_sgRNA_Dynamics)



### > Plot of dynamics ----
# Reshape of the dataset
Reshaped_Simulated_Scenario2_Dynamics <- melt(data=Simulated_Scenario2_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")
Reshaped_Simulated_Real_Dynamics <- melt(data=Simulated_Real_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")

Merged_Simulated_Dynamics_Scenario2 <- rbind(cbind(Type="Counterfactual",Reshaped_Simulated_Scenario2_Dynamics),
                                             cbind(Type="Real",Reshaped_Simulated_Real_Dynamics))
# Calculation of distributions
Distribution_Dynamics_Scenario2 <- Distribution_Function(.data=Merged_Simulated_Dynamics_Scenario2,.variables = .(Type,Group,ObservationType,time))

facet_labels <- c("ACE2/RBD bind. inhibit° \n (AU/mL)","Genomic RNA \n (copies/mL)","Subgenomic RNA \n (copies/mL)","Target cells \n  (T/T0, %)")
names(facet_labels) <- c("Obs_NAb","gRNA_N","sgRNA_N","TN_T0_Percent")

LOD_values <- data.frame(ObservationType=names(facet_labels),Value=c(0.1,LOD_gRNA,LOD_sgRNA,NA))
LOD_values$ObservationType <- factor(LOD_values$ObservationType,levels=names(facet_labels))

Global_Plot_Scenario2 <- Plot_Real_Vs_Counterfactual_Dynamics(data=subset(Distribution_Dynamics_Scenario2,ObservationType %in% c("gRNA_N","sgRNA_N")),facet_labels = facet_labels) + 
  geom_hline(data=subset(LOD_values,ObservationType %in% c("gRNA_N","sgRNA_N")),
             aes(yintercept = Value),color="darkred",linetype="solid",linewidth=0.8,alpha=0.5) +
  scale_x_continuous(name="Time post-exposure (days)",breaks = seq(0,30,by=5))  +
  facetted_pos_scales(y = list(
    ObservationType == "gRNA_N" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(1E-2,NA)),
    ObservationType == "sgRNA_N" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(1E-3,NA)))) +# ,
  coord_cartesian(xlim = c(-0.5,29.5),expand=0) 
# ggsave(Global_Plot_Scenario2,filename=paste(Results_Folder,Figure_Folder,"ViralPlots_Counterfactual_Scenario_Nodepl.png",sep="/"),height=5,width=10,dpi=300)
# ------------ #



# ............................----
# -- COUNTERFACTUAL SCENARIO NO NEUTRALIZATION ----
# ............................----
# Scenario: What if no effect of antibody response on viral dynamics (eta=0) ?


### > Simulation of dynamics ---- 
Simulated_Scenario3_Dynamics <- simulx(model = MlxtranModel_CounterfactualScenario_NoNeut_file,
                                       parameter = Simulated_Population_Parameters,
                                       output = Model_Outputs,regressor = Regressors,
                                       covariate = Covariates,
                                       settings = list(sharedIds="population"))[sapply(Model_Outputs, function(l) l$name)]

Simulated_Scenario3_Dynamics <- Reduce(function(x,y) merge(x,y,all=TRUE), Simulated_Scenario3_Dynamics)
Simulated_Scenario3_Dynamics <- merge(Covariates,Simulated_Scenario3_Dynamics)[-1]
names(Simulated_Scenario3_Dynamics)[which(names(Simulated_Scenario3_Dynamics) == "rep")] <- "pop"
Simulated_Scenario3_Dynamics$Obs_NAb <- 10^Simulated_Scenario3_Dynamics$Obs_NAb
Simulated_Scenario3_Dynamics$gRNA_N <- 10^Simulated_Scenario3_Dynamics$gRNA_N
Simulated_Scenario3_Dynamics$sgRNA_N <- 10^Simulated_Scenario3_Dynamics$sgRNA_N


### > Estimation of VL descriptors ----
Descriptors_Scenario3_gRNA_Dynamics <- ddply(.data=Simulated_Scenario3_Dynamics,.variables = .(pop,Group),summarize,
                                             Peak=Peak_VL_dynamics(time=time,value=gRNA_N),
                                             TimePeak=TimePeak_VL_dynamics(time=time,value=gRNA_N),
                                             AUC=AUC_VL_Dynamics(time=time,value=gRNA_N),
                                             Clearance=Duration_Clearance_Stage(time=time,value=gRNA_N,lod=LOD_gRNA),
                                             Acute=Duration_Acute_Stage(time=time,value=gRNA_N,lod=LOD_gRNA))
Descriptors_Scenario3_gRNA_Dynamics <- cbind(VLType="gRNA",Descriptors_Scenario3_gRNA_Dynamics)

Descriptors_Scenario3_sgRNA_Dynamics <- ddply(.data=Simulated_Scenario3_Dynamics,.variables = .(pop,Group),summarize,
                                              Peak=Peak_VL_dynamics(time=time,value=sgRNA_N),
                                              TimePeak=TimePeak_VL_dynamics(time=time,value=sgRNA_N),
                                              AUC=AUC_VL_Dynamics(time=time,value=sgRNA_N),
                                              Clearance=Duration_Clearance_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA),
                                              Acute=Duration_Acute_Stage(time=time,value=sgRNA_N,lod=LOD_sgRNA))
Descriptors_Scenario3_sgRNA_Dynamics <- cbind(VLType="sgRNA",Descriptors_Scenario3_sgRNA_Dynamics)



### > Plot of dynamics ----
# Reshape of the dataset
Reshaped_Simulated_Scenario3_Dynamics <- melt(data=Simulated_Scenario3_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")
Reshaped_Simulated_Real_Dynamics <- melt(data=Simulated_Real_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")

Merged_Simulated_Dynamics_Scenario3 <- rbind(cbind(Type="Counterfactual",Reshaped_Simulated_Scenario3_Dynamics),
                                             cbind(Type="Real",Reshaped_Simulated_Real_Dynamics))
# Calculation of distributions
Distribution_Dynamics_Scenario3 <- Distribution_Function(.data=Merged_Simulated_Dynamics_Scenario3,.variables = .(Type,Group,ObservationType,time))

facet_labels <- c("ACE2/RBD bind. inhibit° \n (AU/mL)","Genomic RNA \n (copies/mL)","Subgenomic RNA \n (copies/mL)","Target cells \n  (T/T0, %)")
names(facet_labels) <- c("Obs_NAb","gRNA_N","sgRNA_N","TN_T0_Percent")

LOD_values <- data.frame(ObservationType=names(facet_labels),Value=c(0.1,LOD_gRNA,LOD_sgRNA,NA))
LOD_values$ObservationType <- factor(LOD_values$ObservationType,levels=names(facet_labels))

Global_Plot_Scenario3 <- Plot_Real_Vs_Counterfactual_Dynamics(data=subset(Distribution_Dynamics_Scenario3,ObservationType != "Obs_NAb"),facet_labels = facet_labels) + 
  geom_hline(data=subset(LOD_values,ObservationType != "Obs_NAb"),
             aes(yintercept = Value),color="darkred",linetype="solid",linewidth=0.8,alpha=0.5) +
  scale_x_continuous(name="Time post-exposure (days)",breaks = seq(0,30,by=5))  +
  facetted_pos_scales(y = list(
    ObservationType == "gRNA_N" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=4),labels = exp_bold, limits=c(1E-2,5E8)),
    ObservationType == "sgRNA_N" ~ scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=5),labels = exp_bold, limits=c(1E-3,0.5E8)),
    ObservationType == "TN_T0_Percent" ~ scale_y_continuous(breaks = seq(0,100,by=20), limits=c(-1,101)))) +
  coord_cartesian(xlim = c(-0.5,29.5),expand=0) 
# ggsave(Global_Plot_Scenario3,filename=paste(Results_Folder,Figure_Folder,"ViralPlots_Counterfactual_Scenario_NoNeut.png",sep="/"),height=7,width=10,dpi=300)
# ------------ #



# ............................----
# -- MERGED PLOTS SCENARIOS NoDepl & NoNeut ----
# ............................----
library(cowplot)

Legend <- ggpubr::get_legend(Global_Plot_Scenario2)
Merged_plots_Scenario_23_without_Legend <- plot_grid(Global_Plot_Scenario2 + theme(legend.position = "none"),
                                                     Global_Plot_Scenario3 + theme(legend.position = "none"),
                                                     labels = c("A","B"),nrow=2,rel_heights = c(2,3))
Merged_plots_Scenario_23 <- plot_grid(Legend,Merged_plots_Scenario_23_without_Legend,
                                      ncol=1,rel_heights = c(0.05,1))
# ggsave(Merged_plots_Scenario_23,filename=paste(Results_Folder,Figure_Folder,"ViralPlots_Counterfactual_Scenarios_NoNeut_NoDepl.png",sep="/"),height=12,width=10,dpi=300,bg="white")
# ------------ #




# ............................----
# -- PLOT OF VL DESCRIPTORS ----
# ............................----
# Merge of all estimations of descriptors
All_Descriptors <- rbind(cbind(SimulationType="Real",Descriptors_Real_gRNA_Dynamics),
                         cbind(SimulationType="Real",Descriptors_Real_sgRNA_Dynamics),
                         cbind(SimulationType="Scenario1a",Descriptors_Scenario1a_gRNA_Dynamics),
                         cbind(SimulationType="Scenario1a",Descriptors_Scenario1a_sgRNA_Dynamics),
                         cbind(SimulationType="Scenario1b",Descriptors_Scenario1b_gRNA_Dynamics),
                         cbind(SimulationType="Scenario1b",Descriptors_Scenario1b_sgRNA_Dynamics),
                         cbind(SimulationType="Scenario1c",Descriptors_Scenario1c_gRNA_Dynamics),
                         cbind(SimulationType="Scenario1c",Descriptors_Scenario1c_sgRNA_Dynamics),
                         cbind(SimulationType="Scenario2",Descriptors_Scenario2_gRNA_Dynamics),
                         cbind(SimulationType="Scenario2",Descriptors_Scenario2_sgRNA_Dynamics),
                         cbind(SimulationType="Scenario3",Descriptors_Scenario3_gRNA_Dynamics),
                         cbind(SimulationType="Scenario3",Descriptors_Scenario3_sgRNA_Dynamics))

# Calculation of distribution 
Reshaped_All_descriptors <- melt(data=All_Descriptors,id.vars = c("SimulationType","VLType","pop","Group"),variable.name = "Descriptor",value.name = "Value")
Distribution_All_Descriptors <- Distribution_Function(.data=Reshaped_All_descriptors,.variables = .(SimulationType,VLType,Group,Descriptor))


Colors <- Plot_colors_Groups()
Distribution_All_Descriptors$Group <- factor(Distribution_All_Descriptors$Group,levels = names(Colors))
Labels_Descriptors <- setNames(c("Time to peak viral load \n [in days]",
                                 "Peak viral load \n [in copies/mL]",
                                 "Area under viral load curve \n [in copies.day/mL]",
                                 "Duration clearance stage \n [in days]",
                                 "Duration acute stage \n [in days]"), c("TimePeak","Peak","AUC","Clearance","Acute"))
Distribution_All_Descriptors$Descriptor <- factor(Distribution_All_Descriptors$Descriptor,levels=names(Labels_Descriptors))
Distribution_All_Descriptors$SimulationType <- factor(Distribution_All_Descriptors$SimulationType,levels = c("Real","Scenario1a","Scenario1b","Scenario1c","Scenario2","Scenario3"))

strip_settings <- strip_themed(background_y = elem_list_rect(fill=alpha(Colors,alpha = 0.25),color=alpha(Colors,alpha = 1)),by_layer_y = FALSE)

LOD_area <- data.frame(Descriptor=names(Labels_Descriptors),
                       xmax=c(NA,LOD_gRNA,LOD_gRNA*30,NA,NA),xmin=c(NA,0.001,0.001,NA,NA),
                       ymin=c(NA,0,0,NA,NA),ymax=c(NA,7,7,NA,NA))
LOD_area$Descriptor <- factor(LOD_area$Descriptor,levels=names(Labels_Descriptors))


Plot_All_Descriptors_gRNA <- ggplot(data=subset(Distribution_All_Descriptors,VLType == "gRNA")) + 
  geom_rect(data=LOD_area,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=alpha("orange",alpha=0.3)) + 
  geom_vline(data=subset(Distribution_All_Descriptors,SimulationType == "Real" & VLType == "gRNA"),
             aes(xintercept = Median),color="darkred",linetype="dashed",alpha=0.75,linewidth=1.0) +
  geom_pointrange(aes(x=Median,xmin=Q1,xmax=Q3,y=SimulationType,color=Group),shape=21,stroke=1.5,cex=0.60,fill="white",linewidth=1.0,show.legend = FALSE) + 
  geom_pointrange(data=subset(Distribution_All_Descriptors,SimulationType == "Real" & VLType == "gRNA"),
                  aes(x=Median,xmin=Q1,xmax=Q3,y=SimulationType,color=Group),shape=21,stroke=1.5,cex=0.8,fill="darkred",linewidth=1.0,show.legend = FALSE) + 
  
  facet_nested(Group~Descriptor,scales="free",labeller=labeller(Descriptor = Labels_Descriptors),strip=strip_settings) +
  scale_color_manual(values = Colors) + 
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks=seq(0,16,by=2),limits=c(0,16)),
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x,n=6),labels = exp_bold,limits=c(0.5E5,0.5E14)),
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x,n=6),labels = exp_bold,limits=c(0.5E4,1E14)),
    scale_x_continuous(breaks=seq(0,30,by=5),limits=c(0,30)),
    scale_x_continuous(breaks=seq(0,30,by=5),limits=c(0,30))
  )) +
  scale_y_discrete(labels=c("Real","LownAb","MidnAb","HinAb","Nodepl","Noneut")) +
  
  coord_cartesian(ylim=c(1,6)) +
  theme_bw() + 
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(color="black",size=0.75),
        axis.text.x =  element_text(color="black",size=10,face="bold"),
        axis.text.y =  element_text(size=11,face="bold",colour=c("black",rep("gray40",5))),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0,units = "cm"),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(color="black",face = "bold.italic",size = 12),
        strip.background = element_rect(fill=alpha("darkblue",alpha = 0.2),color="darkblue",linewidth = 0.75),
        strip.text.x = element_text(color="black",face="bold",size=11,margin = margin(b=1,t=1)),
        strip.text.y = element_text(color="black",face="bold",size=10.3),
        legend.position = "none") 
# ggsave(Plot_All_Descriptors_gRNA,filename=paste(Results_Folder,Figure_Folder,"MergedPlots_gRNA_Counterfactual_Descriptors.png",sep="/"),height=10,width=12,dpi = 300)





LOD_area <- data.frame(Descriptor=names(Labels_Descriptors),
                       xmax=c(NA,LOD_sgRNA,LOD_sgRNA*30,NA,NA),xmin=c(NA,0.001,0.001,NA,NA),
                       ymin=c(NA,0,0,NA,NA),ymax=c(NA,7,7,NA,NA))
LOD_area$Descriptor <- factor(LOD_area$Descriptor,levels=names(Labels_Descriptors))


Plot_All_Descriptors_sgRNA <- ggplot(data=subset(Distribution_All_Descriptors,VLType == "sgRNA")) + 
  geom_rect(data=LOD_area,aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill=alpha("orange",alpha = 0.3)) + 
  geom_vline(data=subset(Distribution_All_Descriptors,SimulationType == "Real" & VLType == "sgRNA"),
             aes(xintercept = Median),color="darkred",linetype="dashed",alpha=0.75,linewidth=1.0) +
  geom_pointrange(aes(x=Median,xmin=Q1,xmax=Q3,y=SimulationType,color=Group),shape=21,stroke=1.5,cex=0.60,fill="white",linewidth=1.0,show.legend = FALSE) + 
  geom_pointrange(data=subset(Distribution_All_Descriptors,SimulationType == "Real" & VLType == "sgRNA"),
                  aes(x=Median,xmin=Q1,xmax=Q3,y=SimulationType,color=Group),shape=21,stroke=1.5,cex=0.8,fill="darkred",linewidth=1.0,show.legend = FALSE) + 
  
  facet_nested(Group~Descriptor,scales="free",labeller=labeller(Descriptor = Labels_Descriptors),strip=strip_settings) +
  scale_color_manual(values = Colors) + 
  facetted_pos_scales(x = list(
    scale_x_continuous(breaks=seq(0,16,by=2),limits=c(0,15)),
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x,n=6),labels = exp_bold,limits=c(1E-2,NA)),
    scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x,n=6),labels = exp_bold,limits=c(1E-2,NA)),
    scale_x_continuous(breaks=seq(0,30,by=5),limits=c(0,30)),
    scale_x_continuous(breaks=seq(0,30,by=5),limits=c(0,30))
  )) +
  scale_y_discrete(labels=c("Real","LownAb","MidnAb","HinAb","Nodepl","Noneut")) +
  
  coord_cartesian(ylim=c(1,6)) +
  theme_bw() + 
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(color="black",size=0.75),
        axis.text.x =  element_text(color="black",size=10,face="bold"),
        axis.text.y =  element_text(size=11,face="bold",colour=c("black",rep("gray40",5))),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0,units = "cm"),
        axis.title = element_blank(),
        panel.grid.minor = element_blank(),
        plot.title = element_text(color="black",face = "bold.italic",size = 12),
        strip.background = element_rect(fill=alpha("darkblue",alpha = 0.2),color="darkblue",linewidth = 0.75),
        strip.text.x = element_text(color="black",face="bold",size=11,margin = margin(b=1,t=1)),
        strip.text.y = element_text(color="black",face="bold",size=10.3),
        legend.position = "none") 
# ggsave(Plot_All_Descriptors_sgRNA,filename=paste(Results_Folder,Figure_Folder,"MergedPlots_sgRNA_Counterfactual_Descriptors.png",sep="/"),height=10,width=12,dpi = 300)

