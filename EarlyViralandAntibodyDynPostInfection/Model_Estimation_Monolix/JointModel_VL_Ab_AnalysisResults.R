# --------------- # 
# DESCRIPTION: Generation, analysis and visualization of the joint model estimated by Monolix (Generation of results for the paper)
#
# Author: Marie Alexandre 
# R version: 4.2.1
# Monolix Version: 2023R1
# Simulx version: 2023R1

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
Extraction_Parameters_Compartment_FUNCTION <- function(parameters){
  
  names(parameters) <- c(unlist(strsplit(names(parameters),split="_pop",fixed=T)))
  
  Reshaped_parameters <- data.frame(beta = as.numeric(parameters["betas_N"]*1E-5),
                                    delta = as.numeric(parameters["delta_N"]),
                                    P_N = as.numeric(parameters["P_N"]),
                                    P_T = as.numeric(parameters["P_N"]*exp(parameters["fact_P_T"])),
                                    data.frame(as.list(parameters[names(parameters) %notin% c("","betas_N","delta_N","P_N","fact_P_T")])),row.names = NULL)
  
  Groups <- c("Naive","CD40.RBDv - Naive","Conv","CD40.RBDv - Conv","CD40.PanCoV - Conv","mRNA - Conv")
  
  # Addition of group effects
  Parameters_group_naso <- NULL
  Parameters_group_trachea <- NULL
  
  for(group in Groups){
    # group <- Groups
    # - Nasopharynx
    parameters_naso <- data.frame(Group = group,
                                  beta = Reshaped_parameters$beta,
                                  delta = Reshaped_parameters$delta*exp(case_when(group == "Conv" ~ Reshaped_parameters$beta_delta_N_Conv_1,
                                                                                  TRUE ~ 0)),
                                  P = Reshaped_parameters$P_N,
                                  Reshaped_parameters[,c("mu","k","c","cI","alpha_VLSG")],
                                  S0 = Reshaped_parameters$S0*exp(case_when(group == "CD40.RBDv - Naive" ~ Reshaped_parameters$beta_S0_GroupEffect_CD40RBDv,
                                                                            group == "Conv" ~ Reshaped_parameters$beta_S0_GroupEffect_Conv,
                                                                            group == "CD40.RBDv - Conv" ~ Reshaped_parameters$beta_S0_GroupEffect_CD40RBDvConv,
                                                                            group == "CD40.PanCoV - Conv" ~ Reshaped_parameters$beta_S0_GroupEffect_PanCoVConv,
                                                                            group == "mRNA - Conv" ~ Reshaped_parameters$beta_S0_GroupEffect_mRNAConv,
                                                                            TRUE ~ 0)),
                                  gamma = Reshaped_parameters[,c("gamma")],
                                  rho = Reshaped_parameters$rho*exp(case_when(group %in% c("Conv","CD40.RBDv - Conv","CD40.PanCoV - Conv","mRNA - Conv") ~ Reshaped_parameters$beta_rho_Convalescence_1,
                                                                              TRUE ~ 0)),
                                  Reshaped_parameters[,c("Smax","tauS")],
                                  theta_BAb = Reshaped_parameters$theta_BAb*exp(case_when(group == "CD40.RBDv - Naive" ~ Reshaped_parameters$beta_theta_BAb_CD40RBDvNaive_1,
                                                                                          group == "Conv" ~ Reshaped_parameters$beta_theta_BAb_Conv_1,
                                                                                          TRUE ~ 0)),
                                  Reshaped_parameters[,c("alpha_NAb","tauAb")],
                                  eta = Reshaped_parameters$eta*exp(case_when(group == "CD40.RBDv - Conv" ~ Reshaped_parameters$beta_eta_ConvXVacc_CD40,
                                                                              group == "CD40.PanCoV - Conv" ~ Reshaped_parameters$beta_eta_ConvXVacc_PanCov,
                                                                              group == "mRNA - Conv" ~ Reshaped_parameters$beta_eta_ConvXVacc_mRNA,
                                                                              TRUE ~ 0)))
    
    Parameters_group_naso <- rbind(Parameters_group_naso,parameters_naso) 
    
    # - Trachea
    parameters_trachea <- parameters_naso
    parameters_trachea$P <- Reshaped_parameters$P_T
    
    Parameters_group_trachea <- rbind(Parameters_group_trachea,parameters_trachea)
  }
  return(list(Nasopharynx=Parameters_group_naso,Trachea=Parameters_group_trachea))
}

# -- Function of simulations --
Simulation_Individual_Dynamics <- function(project,nrep,output,seed=123456){
  simulation_results <- simulx(project = project,
                               parameter = "mlx_PopUncertainSA", # Simulation of parameters with uncertainty on population parameters 
                               nrep=nrep,
                               output = output,
                               settings = list(sharedIds=c("covariate","regressor"),seed=Seed))
  
  simulation_results <- merge(simulation_results[[output$name]],simulation_results$parameter[,c("rep","id","original_id","Group")])
  
  
  return(simulation_results)
}


# -- Functions "individual dynamics" section --
Individual_Plots_gRNA_sgRNA_FUNCTION <- function(model_data,observed_data,colors,animal_labels){
  
  group <- unique(observed_data$Group)
  
  plot <- ggplot(data=model_data) +
    
    geom_hline(yintercept = 749,linetype="solid",color="darkred",linewidth=1.0,alpha=0.8) +
    
    # Model prediction ---
    geom_ribbon(aes(x=Time,ymin=ICMIN,ymax=ICMAX,fill=as.factor(MeasureType),linetype=as.factor(MeasureType)),alpha=0.80,color="black") +
    geom_ribbon(aes(x=Time,ymin=ICMIN,ymax=ICMAX,linetype=as.factor(MeasureType)),fill=NA,color="black") +
    geom_line(aes(x=Time,y=Value,linetype=as.factor(MeasureType)),linewidth=1.25,color="black") +
    
    scale_fill_manual(name="Model Prediction",breaks=c("gRNA","sgRNA"),values=c("gray30","gray80"),labels=c("Genomic RNA","Subgenomic RNA")) + 
    scale_linetype_manual(name="Model Prediction",breaks=c("gRNA","sgRNA"),values=c("solid","longdash"),labels=c("Genomic RNA","Subgenomic RNA")) + 
    
    # Observed data --- 
    new_scale_fill() + 
    geom_point(data=observed_data,
               aes(x=DayPostExpo,y=CensoredValue,shape=as.factor(MeasureType),fill=as.factor(Censored)),color="black",cex=2.0,stroke=1.0) +  
    scale_shape_manual(name="Observed Data",breaks=c("gRNA","sgRNA"),values=c(21,24),labels=c("Genomic RNA","Subgenomic RNA")) + 
    scale_fill_manual(name="Censored",breaks = c(0,1),values=c(as.character(colors[group]),"white"),labels=c("None","Left-censored")) +
    
    facet_nested(Group ~ SubjectID,labeller = labeller(SubjectID = animal_labels)) +
    
    scale_color_manual(name="Group",breaks = names(colors),values=colors,labels=names(colors),drop=FALSE) + 
    guides(fill = guide_legend(override.aes = list(pch = c(21,21),fill = c("gray40","white"))),
           color = guide_legend(override.aes = list(pch = 21, color="black",fill=colors[group]))) +
    
    scale_x_continuous(breaks = seq(0,30,by=5)) + xlab("Time post-exposure (days)") + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = exp_bold) + ylab("Viral load (copies/mL)") +
    
    coord_cartesian(xlim = c(-0.5,30),ylim=c(1,1.0E10),expand = 0) +
    
    theme_bw() + 
    theme(axis.line = element_line(color="black",linewidth=1.0),
          axis.text.x =  element_text(color="black",size=10,face="bold",hjust = .80),
          axis.text.y =  element_text(color="black",size=10,face="bold"),
          axis.ticks = element_line(color="black",linewidth=1),
          axis.ticks.length = unit(0.1,units = "cm"),
          axis.title = element_text(color="black",size=12,face="bold"),
          legend.key = element_rect(color="white",fill="white",size=0.5),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          legend.key.width = unit(units = "cm",x = 1.5),
          legend.key.height = unit(units = "cm",x = 0.5),
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill=alpha(colors[group],0.50),color=colors[group]),
          strip.text.x = element_text(size = 10,face="bold",margin = margin(b=1,t=0)),
          strip.text.y = element_text(size = 10,face="bold",margin = margin(l=2,r=0)))
  
  return(plot)
}
Individual_Plots_BAb_NAb_FUNCTION <- function(model_data,observed_data,colors,animal_labels){
  
  observed_data <- subset(observed_data,!is.na(CensoredValue))
  observed_data$Censored <- factor(observed_data$Censored,levels = c(0,1))
  
  group <- unique(observed_data$Group)
  
  plot <- ggplot(data=model_data) +
    
    geom_hline(yintercept = 0.1,linetype="solid",color="darkred",linewidth=1.0,alpha=0.8) +
    
    # Model prediction ---
    geom_ribbon(aes(x=Time,ymin=ICMIN,ymax=ICMAX,fill=as.factor(MeasureType),linetype=as.factor(MeasureType)),alpha=0.80,color="black") +
    geom_ribbon(aes(x=Time,ymin=ICMIN,ymax=ICMAX,linetype=as.factor(MeasureType)),fill=NA,color="black") +
    geom_line(aes(x=Time,y=Value,linetype=as.factor(MeasureType)),linewidth=1.25,color="black") +
    
    scale_fill_manual(name="Model Prediction",breaks=c("BAb","NAb"),values=c("gray30","gray80"),labels=c("Binding antibody","Neutralizing antibody")) +
    scale_linetype_manual(name="Model Prediction",breaks=c("BAb","NAb"),values=c("solid","longdash"),labels=c("Binding antibody","Neutralizing antibody")) + 
    
    # Observed data --- 
    new_scale_fill() +
    geom_point(data=observed_data,
               aes(x=DayPostExpo,y=CensoredValue,shape=as.factor(MeasureType),fill=as.factor(Censored)),color="black",cex=2.0,stroke=1.0) +
    scale_shape_manual(name="Observed Data",breaks=c("BAb","NAb"),values=c(21,24),labels=c("Binding antibody","Neutralizing antibody")) +
    scale_fill_manual(name="Censored",breaks = c(0,1),values=c(as.character(colors[group]),"white"),labels=c("None","Left-censored"),drop=FALSE) +
    
    facet_nested(Group ~ SubjectID,labeller=labeller(SubjectID = animal_labels)) +
    
    scale_color_manual(name="Group",breaks = names(colors),values=colors,labels=names(colors),drop=FALSE) +
    guides(fill = guide_legend(override.aes = list(pch = c(21,21),fill = c("gray40","white"))),
           color = guide_legend(override.aes = list(pch = 21, color="black",fill=colors[group]))) +
    
    scale_x_continuous(breaks = seq(0,30,by=5)) + xlab("Time post-exposure (days)") + 
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = exp_bold) + ylab("Viral load (copies/mL)") + 
    
    coord_cartesian(xlim = c(-0.5,30),ylim=c(1e-2,1.0E7),expand = 0) + 
    
    theme_bw() + 
    theme(axis.line = element_line(color="black",size=1.0),
          axis.text =  element_text(color="black",size=10,face="bold"),
          axis.ticks = element_line(color="black",size=1),
          axis.ticks.length = unit(0.1,units = "cm"),
          axis.title = element_text(color="black",size=12,face="bold"),
          legend.key = element_rect(color="white",fill="white",size=0.5),
          legend.position = "right",
          legend.title = element_text(color="black",face = "bold",size=12),
          legend.text = element_text(color="black",face="bold",size=10),
          legend.key.width = unit(units = "cm",x = 1.5),
          legend.key.height = unit(units = "cm",x = 0.5),
          panel.grid.major = element_line(colour = "gray90"),
          panel.grid.minor = element_blank(),
          strip.background = element_rect(fill=alpha(colors[group],0.50),color=colors[group]),
          strip.text.x = element_text(size = 10,face="bold",margin = margin(b=1,t=0)),
          strip.text.y = element_text(size = 10,face="bold",margin = margin(l=2,r=0)))
  
  return(plot)
}

# -- Functions "neutralization quality" section -- 
Epsilon_FUNCTION <- function(eta,NAb){
  
  epsilon <- eta*NAb/(1+eta*NAb)
  return(epsilon)
}

# -- Functions "protection threshold and R0" section -- 
Protective_Threshold_FUNCTION <- function(parameters){
  
  with(parameters,{
    NAb_N <- (beta*T0_N*(mu*P_N-delta)-delta*c)/(c*eta*delta)
    NAb_T <- (beta*T0_T*(mu*P_T-delta)-delta*c)/(c*eta*delta)
    NAb = max(NAb_N,NAb_T)
    return(NAb)
  })
}
Basic_Reproduction_Number_FUNCTION <- function(parameters){
  
  with(parameters,{
    if("NAb_0" %notin% names(parameters)){
      NAb_0 <- alpha_NAb*theta_BAb*S0/(log(2)/tauAb)
    }
    
    R0_N <- beta*T0_N*mu*P_N/( delta*(1+eta*NAb_0)*(c + beta*T0_N/(1+eta*NAb_0) ) ) 
    R0_T <- beta*T0_T*mu*P_T/( delta*(1+eta*NAb_0)*(c + beta*T0_T/(1+eta*NAb_0) ) ) 
    R0 <- max(R0_N,R0_T)
    return(R0)
  })
}
# ------------------ #



# --- LIBRAIRIES ----
library(lixoftConnectors,lib.loc = "C:/Users/marie/AppData/Local/R/win-library/4.4/Monolix2023")

library(RsSimulx)

# Libraries for visual representations
library(ggplot2) ; library(lemon)
library(gridExtra)
library(grid)
library(scales)   # Library for scientific notation in y-axis
library(ggnewscale) 
library(ggpubr)
library(ggh4x)
library(tibble)
library(ggpmisc)
library(cowplot)


# Libraries for parallel calculation
library(parallel) ; library(snow) ; library(doSNOW)

library(reshape2) # Library to reshape dataset
library(plyr)  # Library used to calculate data distribution
library(dplyr)  # Library used for the function case_when
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

# -- Monolix Model -- 
MonolixModel_Folder <- paste(Project_Folder,"Model_Estimation_Monolix",sep="/")
MonolixModel_name <- "Final_JointModel_NHP_DeltaInfection_MONOLIX"
Monolix_Project <- paste(MonolixModel_Folder,paste(MonolixModel_name,"mlxtran",sep="."),sep="/")
Simulx_Project <- paste(MonolixModel_Folder,paste(MonolixModel_name,"_Simulation.smlx",sep=""),sep="/")

# -- Figures ---
Figure_Folder <- "Figures"
# Uncomment to save figures 
# dir.create(paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,sep="/"),recursive = TRUE)

# -- Simulation Code --
MlxtranModel_Simulation_file <- paste(MonolixModel_Folder,"MlxtranCode_NHP_JointModel_VL_Ab_PostInfection_SIMULX.txt",sep="/")

# -- Parameters of simulation ---
Nb_ind <- 250     # Number of individual parameters per NHP for individual dynamics
Nb_pop <- 500     # Number of population parameters for median dynamics
Seed <- 123456    # Seed used to simulated parameters for reproductibility of the work
# ------------ #




# ............................----
# -- VISUAL PREDICTIVE CHECK  ----
# Paper: Figure S8 ----
# ............................----
VPC_Outliers_estimation <- function(percentile_data){
  Splits <- unique(percentile_data$split)
  Outliers_area <- NULL ; Outliers_dots <- NULL
  for(s in 1:length(Splits)){
    data_split <- percentile_data[which(percentile_data$split == Splits[s]),]
    x_axis <- data_split$bins_middles
    
    for(i in 1:length(x_axis)){
      x <- x_axis[i]
      
      # dots
      tmp_results_dots <- NULL
      # median 
      if(between(data_split$empirical_median[i],data_split$theoretical_median_piLower[i],data_split$theoretical_median_piUpper[i])){
        tmp_results_dots <- rbind(tmp_results_dots,data.frame(x=x,y=data_split$empirical_median[i],type="Median",Outlier=0,stringsAsFactors = FALSE))
      }else{
        tmp_results_dots <- rbind(tmp_results_dots,data.frame(x=x,y=data_split$empirical_median[i],type="Median",Outlier=1,stringsAsFactors = FALSE))
      }
      # lower bound
      if(between(data_split$empirical_lower[i],data_split$theoretical_lower_piLower[i],data_split$theoretical_lower_piUpper[i])){
        tmp_results_dots <- rbind(tmp_results_dots,data.frame(x=x,y=data_split$empirical_lower[i],type="Lower",Outlier=0,stringsAsFactors = FALSE))
      }else{
        tmp_results_dots <- rbind(tmp_results_dots,data.frame(x=x,y=data_split$empirical_lower[i],type="Lower",Outlier=1,stringsAsFactors = FALSE))
      }
      # upper bound
      if(between(data_split$empirical_upper[i],data_split$theoretical_upper_piLower[i],data_split$theoretical_upper_piUpper[i])){
        tmp_results_dots <- rbind(tmp_results_dots,data.frame(x=x,y=data_split$empirical_upper[i],type="Upper",Outlier=0,stringsAsFactors = FALSE))
      }else{
        tmp_results_dots <- rbind(tmp_results_dots,data.frame(x=x,y=data_split$empirical_upper[i],type="Upper",Outlier=1,stringsAsFactors = FALSE))
      }
      
      tmp_results_dots$split <- Splits[s]
      Outliers_dots <- rbind(Outliers_dots,tmp_results_dots)
      
      
      # Areas
      if(i < length(x_axis)){
        tmp_results_area <- NULL
        x1 <- x_axis[i] ; x2 <- x_axis[i+1]
        additional_x <- seq(x1,x2,length.out = 500)
        # Estimation of the empirical value for the additional point (linear regression between x1 and x2)
        
        # Empirical median
        Y1_emp_med <- data_split$empirical_median[i] ; Y2_emp_med <- data_split$empirical_median[i+1]  
        Slope_emp_med <- (Y2_emp_med - Y1_emp_med)/(x2-x1) ; Int_emp_med <- Y1_emp_med -  Slope_emp_med*x1
        add_Y_emp_med <- sapply(additional_x,function(new_x,slope,inter){return(slope*new_x+inter)},slope=Slope_emp_med,inter=Int_emp_med)
        # Theoretical median lower bound PI
        Y1_theo_med_piLower <- data_split$theoretical_median_piLower[i] ; Y2_theo_med_piLower <- data_split$theoretical_median_piLower[i+1]
        Slope_theo_med_piLower <- (Y2_theo_med_piLower-Y1_theo_med_piLower)/(x2-x1) ; Int_theo_med_piLower <- Y1_theo_med_piLower - Slope_theo_med_piLower*x1
        add_Y_theo_med_piLower <- sapply(additional_x,function(new_x,slope,inter){return(slope*new_x+inter)},slope=Slope_theo_med_piLower,inter=Int_theo_med_piLower)
        # Theoretical median upper bound PI
        Y1_theo_med_piUpper <- data_split$theoretical_median_piUpper[i] ; Y2_theo_med_piUpper <- data_split$theoretical_median_piUpper[i+1]
        Slope_theo_med_piUpper <- (Y2_theo_med_piUpper-Y1_theo_med_piUpper)/(x2-x1) ; Int_theo_med_piUpper <- Y1_theo_med_piUpper - Slope_theo_med_piUpper*x1
        add_Y_theo_med_piUpper <- sapply(additional_x,function(new_x,slope,inter){return(slope*new_x+inter)},slope=Slope_theo_med_piUpper,inter=Int_theo_med_piUpper)
        
        # Empirical Lower
        Y1_emp_lower <- data_split$empirical_lower[i] ; Y2_emp_lower <- data_split$empirical_lower[i+1]
        Slope_emp_lower <- (Y2_emp_lower - Y1_emp_lower)/(x2-x1) ; Int_emp_lower <- Y1_emp_lower -  Slope_emp_lower*x1
        add_Y_emp_lower <- sapply(additional_x,function(new_x,slope,inter){return(slope*new_x+inter)},slope=Slope_emp_lower,inter=Int_emp_lower)
        # Theoretical Lower lower bound PI
        Y1_theo_lower_piLower <- data_split$theoretical_lower_piLower[i] ; Y2_theo_lower_piLower <- data_split$theoretical_lower_piLower[i+1]
        Slope_theo_lower_piLower <- (Y2_theo_lower_piLower-Y1_theo_lower_piLower)/(x2-x1) ; Int_theo_lower_piLower <- Y1_theo_lower_piLower - Slope_theo_lower_piLower*x1
        add_Y_theo_lower_piLower <- sapply(additional_x,function(new_x,slope,inter){return(slope*new_x+inter)},slope=Slope_theo_lower_piLower,inter=Int_theo_lower_piLower)
        # Theoretical Lower upper bound PI
        Y1_theo_lower_piUpper <- data_split$theoretical_lower_piUpper[i] ; Y2_theo_lower_piUpper <- data_split$theoretical_lower_piUpper[i+1]
        Slope_theo_lower_piUpper <- (Y2_theo_lower_piUpper-Y1_theo_lower_piUpper)/(x2-x1) ; Int_theo_lower_piUpper <- Y1_theo_lower_piUpper - Slope_theo_lower_piUpper*x1
        add_Y_theo_lower_piUpper <- sapply(additional_x,function(new_x,slope,inter){return(slope*new_x+inter)},slope=Slope_theo_lower_piUpper,inter=Int_theo_lower_piUpper)
        
        # Empirical Upper
        Y1_emp_upper <- data_split$empirical_upper[i] ; Y2_emp_upper <- data_split$empirical_upper[i+1]
        Slope_emp_upper <- (Y2_emp_upper - Y1_emp_upper)/(x2-x1) ; Int_emp_upper <- Y1_emp_upper -  Slope_emp_upper*x1
        add_Y_emp_upper <- sapply(additional_x,function(new_x,slope,inter){return(slope*new_x+inter)},slope=Slope_emp_upper,inter=Int_emp_upper)
        # Theoretical upper lower bound PI
        Y1_theo_upper_piLower <- data_split$theoretical_upper_piLower[i] ; Y2_theo_upper_piLower <- data_split$theoretical_upper_piLower[i+1]
        Slope_theo_upper_piLower <- (Y2_theo_upper_piLower-Y1_theo_upper_piLower)/(x2-x1) ; Int_theo_upper_piLower <- Y1_theo_upper_piLower - Slope_theo_upper_piLower*x1
        add_Y_theo_upper_piLower <- sapply(additional_x,function(new_x,slope,inter){return(slope*new_x+inter)},slope=Slope_theo_upper_piLower,inter=Int_theo_upper_piLower)
        # Theoretical Lower upper bound PI
        Y1_theo_upper_piUpper <- data_split$theoretical_upper_piUpper[i] ; Y2_theo_upper_piUpper <- data_split$theoretical_upper_piUpper[i+1]
        Slope_theo_upper_piUpper <- (Y2_theo_upper_piUpper-Y1_theo_upper_piUpper)/(x2-x1) ; Int_theo_upper_piUpper <- Y1_theo_upper_piUpper - Slope_theo_upper_piUpper*x1
        add_Y_theo_upper_piUpper <- sapply(additional_x,function(new_x,slope,inter){return(slope*new_x+inter)},slope=Slope_theo_upper_piUpper,inter=Int_theo_upper_piUpper)
        
        for(j in 1:length(additional_x)){
          # Median
          if(between(add_Y_emp_med[j],add_Y_theo_med_piLower[j],add_Y_theo_med_piUpper[j])){
            tmp_results_area <- rbind(tmp_results_area,data.frame(x=additional_x[j],type="Median",ymin=add_Y_emp_med[j],ymax=add_Y_emp_med[j],Outlier=0,stringsAsFactors = FALSE))
          }else if(add_Y_emp_med[j] > add_Y_theo_med_piUpper[j]){
            tmp_results_area <- rbind(tmp_results_area,data.frame(x=additional_x[j],type="Median",ymin=add_Y_theo_med_piUpper[j],ymax=add_Y_emp_med[j],Outlier=1,stringsAsFactors = FALSE))
          }else{
            tmp_results_area <- rbind(tmp_results_area,data.frame(x=additional_x[j],type="Median",ymin=add_Y_emp_med[j],ymax=add_Y_theo_med_piLower[j],Outlier=1,stringsAsFactors = FALSE))
          }
          
          # Lower
          if(between(add_Y_emp_lower[j],add_Y_theo_lower_piLower[j],add_Y_theo_lower_piUpper[j])){
            tmp_results_area <- rbind(tmp_results_area,data.frame(x=additional_x[j],type="Lower",ymin=add_Y_emp_lower[j],ymax=add_Y_emp_lower[j],Outlier=0,stringsAsFactors = FALSE))
          }else if(add_Y_emp_lower[j] > add_Y_theo_lower_piUpper[j]){
            tmp_results_area <- rbind(tmp_results_area,data.frame(x=additional_x[j],type="Lower",ymin=add_Y_theo_lower_piUpper[j],ymax=add_Y_emp_lower[j],Outlier=1,stringsAsFactors = FALSE))
          }else{
            tmp_results_area <- rbind(tmp_results_area,data.frame(x=additional_x[j],type="Lower",ymin=add_Y_emp_lower[j],ymax=add_Y_theo_lower_piLower[j],Outlier=1,stringsAsFactors = FALSE))
          }
          
          # Upper
          if(between(add_Y_emp_upper[j],add_Y_theo_upper_piLower[j],add_Y_theo_upper_piUpper[j])){
            tmp_results_area <- rbind(tmp_results_area,data.frame(x=additional_x[j],type="Upper",ymin=add_Y_emp_upper[j],ymax=add_Y_emp_upper[j],Outlier=0,stringsAsFactors = FALSE))
          }else if(add_Y_emp_upper[j] > add_Y_theo_upper_piUpper[j]){
            tmp_results_area <- rbind(tmp_results_area,data.frame(x=additional_x[j],type="Upper",ymin=add_Y_theo_upper_piUpper[j],ymax=add_Y_emp_upper[j],Outlier=1,stringsAsFactors = FALSE))
          }else{
            tmp_results_area <- rbind(tmp_results_area,data.frame(x=additional_x[j],type="Upper",ymin=add_Y_emp_upper[j],ymax=add_Y_theo_upper_piLower[j],Outlier=1,stringsAsFactors = FALSE))
          }
        }
        tmp_results_area$split <- Splits[s]
        Outliers_area <- rbind(Outliers_area,tmp_results_area)
      }
    }
  }
  return(list(Outliers_dots=Outliers_dots,Outliers_area=Outliers_area))
}
VPC_GroupName_change <- function(data,monolix_names){
  data$split <- revalue(data$split,setNames(names(Plot_colors_Groups()),paste("GroupEffect: ",Monolix_Group_names,sep="")))
  return(data)
}
VPC_log10_Scale <- function(lab) {
  do.call(
    expression,
    lapply(paste(lab), function(x) bquote(bold("10"^.(x))))
  )
}


# Load of the project
initializeLixoftConnectors(software = "monolix",path = "C:/ProgramData/Lixoft/MonolixSuite2023R1",force=TRUE)
loadProject(projectFile = Monolix_Project)

Monolix_Group_names <- c("Naive","CD40RBDv","Conv","CD40RBDvConv","PanCoVConv","mRNAConv")
Colors <- Plot_colors_Groups()

# Load of data
Outliers_dots_data <- NULL
Outliers_area_data <- NULL

## > Genomic RNA in nasopharynx ----
VPC_data_gRNAN <- getChartsData(plotName = "plotVpc",
                                computeSettings = list(obsName="y1",
                                                       xBinsSettings = list(is.fixedNbBins=TRUE,nbBins=10)),
                                splitGroup = list(name="GroupEffect"))
# Modification of group names 
VPC_data_gRNAN <- lapply(VPC_data_gRNAN, function(list) VPC_GroupName_change(data=list,monolix_names = Monolix_Group_names))
# Estimation of the outliers (information not included in Monolix files)
Outliers_data_gRNAN <- VPC_Outliers_estimation(percentile_data = VPC_data_gRNAN$percentiles)
# ---

## > Subgenomic RNA in nasopharynx ----
VPC_data_sgRNAN <- getChartsData(plotName = "plotVpc",
                                 computeSettings = list(obsName="y2",
                                                        xBinsSettings = list(is.fixedNbBins=TRUE,nbBins=10)),
                                 splitGroup = list(name="GroupEffect"))
# Modification of group names 
VPC_data_sgRNAN <- lapply(VPC_data_sgRNAN, function(list) VPC_GroupName_change(data=list,monolix_names = Monolix_Group_names))
# Estimation of the outliers (information not included in Monolix files)
Outliers_data_sgRNAN <- VPC_Outliers_estimation(percentile_data = VPC_data_sgRNAN$percentiles)
# ---

## > Genomic RNA in trachea ----
VPC_data_gRNAT <- getChartsData(plotName = "plotVpc",
                                computeSettings = list(obsName="y3",
                                                       xBinsSettings = list(is.fixedNbBins=TRUE,nbBins=10)),
                                splitGroup = list(name="GroupEffect"))
# Modification of group names 
VPC_data_gRNAT <- lapply(VPC_data_gRNAT, function(list) VPC_GroupName_change(data=list,monolix_names = Monolix_Group_names))
# Estimation of the outliers (information not included in Monolix files)
Outliers_data_gRNAT <- VPC_Outliers_estimation(percentile_data = VPC_data_gRNAT$percentiles)
# ---

## > Subgenomic RNA in trachea ----
VPC_data_sgRNAT <- getChartsData(plotName = "plotVpc",
                                 computeSettings = list(obsName="y4",
                                                        xBinsSettings = list(is.fixedNbBins=TRUE,nbBins=10)),
                                 splitGroup = list(name="GroupEffect"))
# Modification of group names 
VPC_data_sgRNAT <- lapply(VPC_data_sgRNAT, function(list) VPC_GroupName_change(data=list,monolix_names = Monolix_Group_names))
# Estimation of the outliers (information not included in Monolix files)
Outliers_data_sgRNAT <- VPC_Outliers_estimation(percentile_data = VPC_data_sgRNAT$percentiles)
# ---

## > Binding antibodies ----
VPC_data_BAb <- getChartsData(plotName = "plotVpc",
                              computeSettings = list(obsName="y5",
                                                     xBinsSettings = list(is.fixedNbBins=TRUE,nbBins=5)),
                              splitGroup = list(name="GroupEffect"))
# Modification of group names 
VPC_data_BAb <- lapply(VPC_data_BAb, function(list) VPC_GroupName_change(data=list,monolix_names = Monolix_Group_names))
# Estimation of the outliers (information not included in Monolix files)
Outliers_data_BAb <- VPC_Outliers_estimation(percentile_data = VPC_data_BAb$percentiles)
# ---

## > Neutralizing antibodies ----
VPC_data_NAb <- getChartsData(plotName = "plotVpc",
                              computeSettings = list(obsName="y6",
                                                     xBinsSettings = list(is.fixedNbBins=TRUE,nbBins=5)),
                              splitGroup = list(name="GroupEffect"))
# Modification of group names 
VPC_data_NAb <- lapply(VPC_data_NAb, function(list) VPC_GroupName_change(data=list,monolix_names = Monolix_Group_names))
# Estimation of the outliers (information not included in Monolix files)
Outliers_data_NAb <- VPC_Outliers_estimation(percentile_data = VPC_data_NAb$percentiles)
# ---

# Merge of data
Percentile_Data <- rbind(cbind(obsName = "gRNAN",VPC_data_gRNAN$percentiles),
                         cbind(obsName = "sgRNAN",VPC_data_sgRNAN$percentiles),
                         cbind(obsName = "gRNAT",VPC_data_gRNAT$percentiles),
                         cbind(obsName = "sgRNAT",VPC_data_sgRNAT$percentiles),
                         cbind(obsName = "BAb",VPC_data_BAb$percentiles),
                         cbind(obsName = "NAb",VPC_data_NAb$percentiles))
Percentile_Data$obsName <- factor(Percentile_Data$obsName,levels=c("gRNAN","sgRNAN","gRNAT","sgRNAT","BAb","NAb"))

Outliers_dots_data <- rbind(cbind(obsName = "gRNAN",Outliers_data_gRNAN$Outliers_dots) ,
                            cbind(obsName = "sgRNAN",Outliers_data_sgRNAN$Outliers_dots),
                            cbind(obsName = "gRNAT",Outliers_data_gRNAT$Outliers_dots),
                            cbind(obsName = "sgRNAT",Outliers_data_sgRNAT$Outliers_dots),
                            cbind(obsName = "BAb",Outliers_data_BAb$Outliers_dots),
                            cbind(obsName = "NAb",Outliers_data_NAb$Outliers_dots))
Outliers_dots_data$obsName <- factor(Outliers_dots_data$obsName,levels=c("gRNAN","sgRNAN","gRNAT","sgRNAT","BAb","NAb"))

Outliers_area_data <- rbind(cbind(obsName = "gRNAN",Outliers_data_gRNAN$Outliers_area) ,
                            cbind(obsName = "sgRNAN",Outliers_data_sgRNAN$Outliers_area),
                            cbind(obsName = "gRNAT",Outliers_data_gRNAT$Outliers_area),
                            cbind(obsName = "sgRNAT",Outliers_data_sgRNAT$Outliers_area),
                            cbind(obsName = "BAb",Outliers_data_BAb$Outliers_area),
                            cbind(obsName = "NAb",Outliers_data_NAb$Outliers_area))
Outliers_area_data <- Outliers_area_data[!duplicated(Outliers_area_data),]
Outliers_area_data$obsName <- factor(Outliers_area_data$obsName,levels=c("gRNAN","sgRNAN","gRNAT","sgRNAT","BAb","NAb"))




# VPC plots 
strip_settings <- strip_themed(background_x = elem_list_rect(fill=alpha(Colors,alpha = 0.25),color=alpha(Colors,alpha = 1)),
                               by_layer_x = FALSE,
                               background_y = element_rect(fill="white",color="white"),
                               text_x = element_text(color="black",size=10,face="bold"),
                               text_y = element_text(color="black",size=11,face="bold"))

Facet_y_labels <- setNames(c("gRNA Naso. \n (copies/mL)","sgRNA Naso. \n (copies/mL)",
                             "gRNA Trachea \n (copies/mL)","sgRNA Trachea \n (copies/mL)",
                             "Anti-RBD IgG \n (AU/mL)","ACE2/RBD bind. \n inhibit° (AU/mL)"),
                           c("gRNAN","sgRNAN","gRNAT","sgRNAT","BAb","NAb"))

VPC <- ggplot(data=Percentile_Data,aes(x=bins_middles)) + 
  # Outliers areas
  geom_ribbon(data=Outliers_area_data,aes(x=x,ymin=ymin,ymax=ymax,group=as.factor(type)),alpha=0.8,fill="#ff0000") + 
  
  # Lower percentile
  geom_ribbon(aes(ymin=theoretical_lower_piLower,ymax=theoretical_lower_piUpper),alpha=0.15,fill=rgb(0/250,128/250,255/255,1)) + 
  geom_line(aes(y=empirical_lower),size=1.0,color="darkblue",linetype="dashed") +
  geom_point(data=subset(Outliers_dots_data,Outlier==0 & type == "Lower"),aes(y=y,x=x),shape=21,fill="white",cex=1.5) + 
  
  # Upper percentile
  geom_ribbon(aes(ymin=theoretical_upper_piLower,ymax=theoretical_upper_piUpper),alpha=0.15,fill=rgb(0/250,128/250,255/255,1)) + 
  geom_line(aes(y=empirical_upper),size=1.0,color="darkblue",linetype="dashed") +
  geom_point(data=subset(Outliers_dots_data,Outlier==0 & type == "Upper"),aes(y=y,x=x),shape=21,fill="white",cex=1.5) +
  
  # Median
  geom_ribbon(aes(ymin=theoretical_median_piLower,ymax=theoretical_median_piUpper),alpha=0.4,fill="navyblue") +
  geom_line(aes(y=empirical_median),size=1.0,color="black") +
  geom_point(data=subset(Outliers_dots_data,Outlier==0 & type == "Median"),aes(y=y,x=x),shape=21,fill="white",cex=1.5) +
  
  # Outliers dots
  geom_point(data=subset(Outliers_dots_data,Outlier == 1),aes(x=x,y=y),fill="red",color="black",pch=21,stroke=1.05,cex=1.5,show.legend = FALSE) +
  
  facet_nested(obsName~split,scales = "free_y",switch="y",strip=strip_settings,labeller=labeller(obsName=Facet_y_labels)) +
  
  theme_bw() +
  theme(
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(colour = "gray98"),
    axis.line = element_line(color="black",size=1),
    axis.text =  element_text(color="black",size=9,face="bold"),
    axis.title = element_text(color="black",size=11,face="bold"),
    panel.spacing.x = unit(0,units="cm")) +
  theme(strip.placement.y = "outside",
        strip.text.y = element_text(margin = margin(unit = "cm",r = 0,l=0)),
        axis.title.y = element_blank()) + 
  
  xlab("Time post-exposure (days)") + scale_x_continuous(breaks = seq(0,29,by=5)) + 
  facetted_pos_scales(y = list(
    obsName == "gRNAN" ~ scale_y_continuous(breaks = seq(-6,10,by=2),labels = VPC_log10_Scale,limits = c(-6,10)),
    obsName == "sgRNAN" ~ scale_y_continuous(breaks = seq(-6,8,by=2),labels = VPC_log10_Scale,limits = c(-6,8)),
    obsName == "gRNAT" ~ scale_y_continuous(breaks = seq(-6,10,by=2),labels = VPC_log10_Scale,limits = c(-6,10)),
    obsName == "sgRNAT" ~ scale_y_continuous(breaks = seq(-6,8,by=2),labels = VPC_log10_Scale,limits = c(-6,8)),
    obsName == "BAb" ~ scale_y_continuous(breaks = seq(1,6,by=1),labels = VPC_log10_Scale,limits = c(1,6)),
    obsName == "NAb" ~ scale_y_continuous(breaks = seq(-2,3,by=1),labels = VPC_log10_Scale,limits = c(-2.5,3))
  ))
# ggsave(VPC,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"Visual_Predictive_Check.png",sep="/"),height=10,width=10,dpi = 300)







# ............................----
# -- CONVERGENCE ASSESSMENT (visualization) ----
# Paper: Figures S9 and S10 ----
# ............................----
# /!\/!\ Convergence assessment should have been run before running the following code to visualize results  /!\/!\
LowerBound_CI <- function(meanValue,sdValue,distribution="LN"){
  
  if(distribution == "LN"){
    icmin <- exp(log(meanValue) - 1.036*sdValue/meanValue)
  }else if(distribution == "N"){
    icmin <- meanValue - 1.036*sdValue
  }else{
    icmin <- NA
  }
  return(icmin)
}
UpperBound_CI <- function(meanValue,sdValue,distribution="LN"){
  
  if(distribution == "LN"){
    icmin <- exp(log(meanValue) + 1.036*sdValue/meanValue)
  }else if(distribution == "N"){
    icmin <- meanValue + 1.036*sdValue
  }else{
    icmin <- NA
  }
  return(icmin)
}


# Load of the project
initializeLixoftConnectors(software = "monolix",path = "C:/ProgramData/Lixoft/MonolixSuite2023R1",force=TRUE)
loadProject(projectFile = Monolix_Project)


## Study of parameter distributions  ----
### > Extraction of results ----
Nb_Run <- 10  # To modify according to the number of iterations ran.
SubProject_Names <- paste(MonolixModel_Folder,MonolixModel_name,"Assessment",
                          paste(MonolixModel_name,stringr::str_pad(seq(1,Nb_Run),width=3,side="left",pad ="0"),sep="_"),
                          paste(MonolixModel_name,"mlxtran",sep="."),sep="/")



Population_Parmeters_Estimation <- NULL
Loglikelihood <- NULL
for(i in 1:Nb_Run){
  # i <- 1
  subproject_name <- SubProject_Names[i]
  
  # Load project
  loadProject(subproject_name)
  # Extraction of population parameters
  estimated_population_parameters <- getEstimatedPopulationParameters()
  population_parameters <- data.frame(parameter=names(estimated_population_parameters),Value=estimated_population_parameters,row.names = NULL)
  estimated_standardError <- getEstimatedStandardErrors()$stochasticApproximation
  estimated_standardError <- estimated_standardError[,-ncol(estimated_standardError)]
  population_parameters <- merge(population_parameters,estimated_standardError,by="parameter")
  population_parameters <- population_parameters[match(estimated_standardError$parameter,population_parameters$parameter),]
  population_parameters$rse <- abs(population_parameters$se*100/population_parameters$Value)
  
  # Extraction of LL
  loglikelihood <- getEstimatedLogLikelihood()$importanceSampling
  Loglikelihood <- rbind(Loglikelihood,data.frame(Run=i,parameter="OFV",Value=loglikelihood["OFV"],se=NA,rse=NA))
  
  Population_Parmeters_Estimation <- rbind(Population_Parmeters_Estimation,cbind(Run=i,population_parameters))
}
# Estimation of confidence interval
Population_Parmeters_Estimation$ICMIN <- with(Population_Parmeters_Estimation,
                                              case_when(parameter %in% c("betas_N_pop","delta_N_pop","P_N_pop","alpha_VLSG_pop","rho_pop","Smax_pop","theta_BAb_pop","alpha_NAb_pop","eta_pop") ~ LowerBound_CI(meanValue = Value,sdValue = se,distribution="LN"),
                                                        TRUE ~ LowerBound_CI(meanValue = Value,sdValue = se,distribution="N")))
Population_Parmeters_Estimation$ICMAX <- with(Population_Parmeters_Estimation,
                                              case_when(parameter %in% c("betas_N_pop","delta_N_pop","P_N_pop","alpha_VLSG_pop","rho_pop","Smax_pop","theta_BAb_pop","alpha_NAb_pop","eta_pop") ~ UpperBound_CI(meanValue = Value,sdValue = se,distribution="LN"),
                                                        TRUE ~ UpperBound_CI(meanValue = Value,sdValue = se,distribution="N")))
Population_Parmeters_Estimation <- rbind(Population_Parmeters_Estimation,cbind(Loglikelihood,ICMIN=NA,ICMAX=NA))
Population_Parmeters_Estimation$parameter <- factor(Population_Parmeters_Estimation$parameter,levels=unique(Population_Parmeters_Estimation$parameter))


### > Plot of results ----
# Modify the order according to parameters ranking
Labels_Pop_params <- c("beta0 (x10^5)","delta","phi_delta_Conv","PN","f_T_P","alpha_sgRNA",
                       "phi_S0_CD40-Naive","Phi_S0_CD40-Conv","phi_S0_mRNA-Conv","phi_S0_PanCoV-Conv",
                       "rho","phi_rho_Conval","Smax","theta","phi_theta_CD40-Naive","phi_theta_Conv",
                       "alpha","eta","phi_eta_CD40-Conv","phi_eta_PanCoV-Conv","phi_eta_mRNA-Conv",
                       "omega_delta","omega_S0","omega_theta","omega_eta",
                       "sigma1T","sigma2T","sigma1N","sigma2N","sigma3","sigma4","-2LL")
names(Labels_Pop_params) <- unique(Population_Parmeters_Estimation$parameter)


Plot1_PopParameters_Distribution <- ggplot(data=subset(Population_Parmeters_Estimation,parameter %in% names(Labels_Pop_params)[1:18])) + 
  # lower limit 
  geom_segment(aes(x=Run-0.30,xend=Run+0.30,y=ICMIN,yend=ICMIN,color=as.factor(Run)),linewidth = 0.80) + 
  # upper limit
  geom_segment(aes(x=Run-0.30,xend=Run+0.30,y=ICMAX,yend=ICMAX,color=as.factor(Run)),linewidth = 0.80) + 
  # Value
  geom_point(aes(x=Run,y=Value,color=as.factor(Run)),cex=1.5,pch=19) + 
  
  facet_rep_wrap(.~parameter,repeat.tick.labels = TRUE,labeller=labeller(parameter=Labels_Pop_params),scales = "free_y",ncol=3) +
  scale_x_continuous(breaks=seq(1,Nb_Run)) + 
  
  theme_bw() + 
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.05,units = "cm"),
        axis.title = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill=alpha("darkblue",alpha = 0.2),color="darkblue",linewidth = 0.75),
        strip.text = element_text(color="black",face="bold",size=10,margin = margin(b=2,t=0)))
# ggsave(Plot1_PopParameters_Distribution,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"ConvergenceAssessment_ParameterDistributions1.png",sep="/"),dpi = 300,height=24,width=18,units = "cm")

Plot2_PopParameters_Distribution <- ggplot(data=subset(Population_Parmeters_Estimation,parameter %in% names(Labels_Pop_params)[19:32])) + 
  # lower limit 
  geom_segment(aes(x=Run-0.30,xend=Run+0.30,y=ICMIN,yend=ICMIN,color=as.factor(Run)),linewidth = 0.80) + 
  # upper limit
  geom_segment(aes(x=Run-0.30,xend=Run+0.30,y=ICMAX,yend=ICMAX,color=as.factor(Run)),linewidth = 0.80) + 
  # Value
  geom_point(aes(x=Run,y=Value,color=as.factor(Run)),cex=1.5,pch=19) + 
  
  facet_rep_wrap(.~parameter,repeat.tick.labels = TRUE,labeller=labeller(parameter=Labels_Pop_params),scales = "free_y",ncol=3) +
  scale_x_continuous(breaks=seq(1,Nb_Run)) + 
  
  theme_bw() + 
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.05,units = "cm"),
        axis.title = element_blank(),
        legend.position = "none",
        panel.grid.minor = element_blank(),
        strip.background = element_rect(fill=alpha("darkblue",alpha = 0.2),color="darkblue",linewidth = 0.75),
        strip.text = element_text(color="black",face="bold",size=10,margin = margin(b=2,t=0)))
# ggsave(Plot2_PopParameters_Distribution,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"ConvergenceAssessment_ParameterDistributions2.png",sep="/"),height=20,width=18,unit="cm",dpi = 300)





# ........................----
# -- INDIVIDUAL DYNAMICS  ----
# Paper: Figures S6, S7 ----
# ........................----
# Load of the project
initializeLixoftConnectors(software = "monolix",path = "C:/ProgramData/Lixoft/MonolixSuite2023R1",force=TRUE)
loadProject(projectFile = Monolix_Project)

Estimated_Individual_Parameters <- getEstimatedIndividualParameters()$conditionalMode
Time_simulation <- seq(0,30,0.25)
Model_Outputs <- list(list(name="gRNA_N",time=Time_simulation),
                      list(name="gRNA_T",time=Time_simulation),
                      list(name="sgRNA_N",time=Time_simulation),
                      list(name="sgRNA_T",time=Time_simulation),
                      list(name="Obs_BAb",time=Time_simulation),
                      list(name="Obs_NAb",time=Time_simulation))

initializeLixoftConnectors(software = "simulx",path = "C:/ProgramData/Lixoft/MonolixSuite2023R1",force=TRUE)

Individual_Regressors <- unique(ObservedData[,c("SubjectID","Inoc","Weight")])
Individual_Regressors <- merge(expand.grid(SubjectID=unique(ObservedData$SubjectID),Time=Time_simulation),Individual_Regressors)
Individual_Regressors <- Individual_Regressors[order(Individual_Regressors$Time,Individual_Regressors$SubjectID),]
names(Individual_Regressors) <- c("id","time","Inoc","Weight")

Estimated_Individual_Dynamics <- simulx(project=Monolix_Project,parameter = Estimated_Individual_Parameters,
                                        regressor = Individual_Regressors,output = Model_Outputs)
Estimated_Individual_Dynamics <- Reduce(function(x,y) merge(x, y, all=TRUE,by=c("id","original_id","time")),Estimated_Individual_Dynamics[1:6])
Estimated_Individual_Dynamics <- melt(Estimated_Individual_Dynamics,id.vars = c("id","original_id","time"),variable.name = "ObservationType",value.name = "log10Value")
Estimated_Individual_Dynamics <- Estimated_Individual_Dynamics[,-1]
names(Estimated_Individual_Dynamics) <- c("SubjectID","Time","ObservationType","log10Value")
Estimated_Individual_Dynamics <- merge(unique(ObservedData[,c("SubjectID","Group")]),Estimated_Individual_Dynamics)
Estimated_Individual_Dynamics$Value <- 10^Estimated_Individual_Dynamics$log10Value
Estimated_Individual_Dynamics$MeasureType <- with(Estimated_Individual_Dynamics,
                                                  case_when(ObservationType %in% c("gRNA_N","gRNA_T") ~ "gRNA",
                                                            ObservationType %in% c("sgRNA_N","sgRNA_T") ~ "sgRNA",
                                                            ObservationType == "Obs_BAb" ~ "BAb",
                                                            ObservationType == "Obs_NAb" ~ "NAb"))




## Simulation of individual dynamics (accounting for population uncertainty) ----
initializeLixoftConnectors(software = "simulx",path = "C:/ProgramData/Lixoft/MonolixSuite2023R1",force=TRUE)

Time_simulation <- seq(0,30,0.25)
Model_Outputs <- list(list(name="gRNA_N",time=Time_simulation),
                      list(name="gRNA_T",time=Time_simulation),
                      list(name="sgRNA_N",time=Time_simulation),
                      list(name="sgRNA_T",time=Time_simulation),
                      list(name="Obs_BAb",time=Time_simulation),
                      list(name="Obs_NAb",time=Time_simulation))

# Simulation of parameters and dynamics directly from Monolix Project (inclusion Uncertainty on Population)
{
  cat("Simulation gRNA N \n")
  Simulation_gRNA_N <-  Simulation_Individual_Dynamics(project = Monolix_Project,
                                                       nrep=Nb_ind, seed = Seed,
                                                       output = Model_Outputs[[1]])
  
  cat("Simulation gRNA T \n")
  Simulation_gRNA_T <- Simulation_Individual_Dynamics(project = Monolix_Project,
                                                      nrep=Nb_ind, seed = Seed,
                                                      output = Model_Outputs[[2]])
  
  cat("Simulation sgRNA N \n")
  Simulation_sgRNA_N <-  Simulation_Individual_Dynamics(project = Monolix_Project,
                                                        nrep=Nb_ind, seed = Seed,
                                                        output = Model_Outputs[[3]])  
  cat("Simulation sgRNA T \n")
  Simulation_sgRNA_T <-  Simulation_Individual_Dynamics(project = Monolix_Project,
                                                        nrep=Nb_ind, seed = Seed,
                                                        output = Model_Outputs[[4]])  
  cat("Simulation Obs BAb \n")
  Simulation_Obs_BAb <-  Simulation_Individual_Dynamics(project = Monolix_Project,
                                                        nrep=Nb_ind, seed = Seed,
                                                        output = Model_Outputs[[5]])  
  cat("Simulation Obs NAb \n")
  Simulation_Obs_NAb <-  Simulation_Individual_Dynamics(project = Monolix_Project,
                                                        nrep=Nb_ind, seed = Seed,
                                                        output = Model_Outputs[[6]])
}
Simulated_Individual_dynamics <- Reduce(function(x,y) merge(x, y, all=TRUE,by=c("rep","id","original_id","Group","time")), 
                                        list(Simulation_gRNA_N,Simulation_gRNA_T,Simulation_sgRNA_N,
                                             Simulation_sgRNA_T,Simulation_Obs_BAb,Simulation_Obs_NAb))

Simulated_Individual_dynamics <- melt(Simulated_Individual_dynamics,id.vars = c("rep","id","original_id","Group","time"),variable.name = "ObservationType",value.name = "log10Value")
Simulated_Individual_dynamics$Value <- 10^Simulated_Individual_dynamics$log10Value
names(Simulated_Individual_dynamics)[names(Simulated_Individual_dynamics) == "original_id"] <- "SubjectID"
names(Simulated_Individual_dynamics)[names(Simulated_Individual_dynamics) == "time"]  <- "Time"



## Plot of dynamics ----
# Estimation of individual dynamics distribution (Mean, CI95%)
Distribution_Individual_dynamics <- ddply(.data=Simulated_Individual_dynamics,.variables = .(Group,id,SubjectID,Time,ObservationType),summarize,
                                          Mean = mean(Value,na.rm=TRUE),
                                          Median = median(Value,na.rm=TRUE),
                                          ICMIN = quantile(Value,na.rm=TRUE,probs=c(0.025)),
                                          ICMAX = quantile(Value,na.rm=TRUE,probs=c(0.975)))

Distribution_Individual_dynamics$MeasureType <- with(Distribution_Individual_dynamics,
                                                     case_when(ObservationType %in% c("gRNA_N","gRNA_T") ~ "gRNA",
                                                               ObservationType %in% c("sgRNA_N","sgRNA_T") ~ "sgRNA",
                                                               ObservationType == "Obs_BAb" ~ "BAb",
                                                               ObservationType == "Obs_NAb" ~ "NAb"))

Distribution_Individual_dynamics <- merge(Distribution_Individual_dynamics,Estimated_Individual_Dynamics)

Colors <- Plot_colors_Groups()
Colors["Naive"] <- "#000000"  # Can be directly modified in the function
# Selection of 3 NHPs per groups (for paper figure)
Selected_NHPs <- unlist(tapply(NHPs$SubjectID,NHPs$Group,sample,size=3),use.names = FALSE)
Labels_NPHs <- Selected_NHPs


### > Viral load in Nasopharynx ---- 
# --- Naive group
Plot_VL_Naso_Naive <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "Naive" & ObservationType %in% c("gRNA_N") & SubjectID %in% Selected_NHPs),
                                                           observed_data=subset(ObservedData,Group == "Naive" & MeasureType %in% c("gRNA") & SampleType == "Nasal fluid" & SubjectID %in% Selected_NHPs),
                                                           colors = Colors,animal_labels=Labels_NPHs) 
# --- CD40.RBDv-Naive Group
Plot_VL_Naso_CD40RBDvNaive <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "CD40.RBDv - Naive" & ObservationType %in% c("gRNA_N") & SubjectID %in% Selected_NHPs),
                                                                   observed_data=subset(ObservedData,Group == "CD40.RBDv - Naive" & MeasureType %in% c("gRNA") & SampleType == "Nasal fluid" & SubjectID %in% Selected_NHPs),
                                                                   colors = Colors,animal_labels=Labels_NPHs) 

# --- Conv Group
Plot_VL_Naso_Conv <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "Conv" & ObservationType %in% c("gRNA_N") & SubjectID %in% Selected_NHPs),
                                                          observed_data=subset(ObservedData,Group == "Conv" & MeasureType %in% c("gRNA") & SampleType == "Nasal fluid" & SubjectID %in% Selected_NHPs),
                                                          colors = Colors,animal_labels=Labels_NPHs)

# --- CD40.RBDv-Conv group
Plot_VL_Naso_CD40RBDvConv <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "CD40.RBDv - Conv" & ObservationType %in% c("gRNA_N") & SubjectID %in% Selected_NHPs),
                                                                  observed_data=subset(ObservedData,Group == "CD40.RBDv - Conv" & MeasureType %in% c("gRNA") & SampleType == "Nasal fluid" & SubjectID %in% Selected_NHPs),
                                                                  colors = Colors,animal_labels=Labels_NPHs)

# --- CD40.PanCov-Conv group
Plot_VL_Naso_CD40PanCovConv <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "CD40.PanCoV - Conv" & ObservationType %in% c("gRNA_N") & SubjectID %in% Selected_NHPs),
                                                                    observed_data=subset(ObservedData,Group == "CD40.PanCoV - Conv" & MeasureType %in% c("gRNA") & SampleType == "Nasal fluid" & SubjectID %in% Selected_NHPs),
                                                                    colors = Colors,animal_labels=Labels_NPHs)

# --- mRNA-Conv group
Plot_VL_Naso_mRNAConv <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "mRNA - Conv" & ObservationType %in% c("gRNA_N") & SubjectID %in% Selected_NHPs),
                                                              observed_data=subset(ObservedData,Group == "mRNA - Conv" & MeasureType %in% c("gRNA") & SampleType == "Nasal fluid" & SubjectID %in% Selected_NHPs),
                                                              colors = Colors,animal_labels=Labels_NPHs)


# --- Merged groups 
Legend <- cowplot::get_legend(Plot_VL_Naso_Naive)

Figure <- grid.arrange(Plot_VL_Naso_Naive + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_VL_Naso_CD40RBDvNaive + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_VL_Naso_Conv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_VL_Naso_CD40RBDvConv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_VL_Naso_CD40PanCovConv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_VL_Naso_mRNAConv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       layout_matrix = rbind(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3)))


Plot_VL_Naso_AllGroups <- ggarrange(annotate_figure(Figure, left = textGrob("Viral load (copies/mL)", rot = 90, just="centre", gp = gpar(fontsize = 12,fontface="bold")),
                                                    bottom = textGrob("Time post-exposure (days)",just="centre", gp = gpar(fontsize = 12,fontface="bold")))) + 
  bgcolor("white") + border("white")  
# ggsave(Plot_VL_Naso_AllGroups,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"IndivFits_VL_Naso.png",sep="/"),height=12,width=10,dpi=300)




### > Viral load in Trachea ---- 
# --- Naive group
Plot_VL_Trachea_Naive <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "Naive" & ObservationType %in% c("gRNA_T") & SubjectID %in% Selected_NHPs),
                                                              observed_data=subset(ObservedData,Group == "Naive" & MeasureType %in% c("gRNA") & SampleType == "Tracheal fluid" & SubjectID %in% Selected_NHPs),
                                                              colors = Colors,animal_labels=Labels_NPHs)

# --- CD40.RBDv-Naive Group
Plot_VL_Trachea_CD40RBDvNaive <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "CD40.RBDv - Naive" & ObservationType %in% c("gRNA_T") & SubjectID %in% Selected_NHPs),
                                                                      observed_data=subset(ObservedData,Group == "CD40.RBDv - Naive" & MeasureType %in% c("gRNA") & SampleType == "Tracheal fluid" & SubjectID %in% Selected_NHPs),
                                                                      colors = Colors,animal_labels=Labels_NPHs)

# --- Conv Group
Plot_VL_Trachea_Conv <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "Conv" & ObservationType %in% c("gRNA_T") & SubjectID %in% Selected_NHPs),
                                                             observed_data=subset(ObservedData,Group == "Conv" & MeasureType %in% c("gRNA") & SampleType == "Tracheal fluid" & SubjectID %in% Selected_NHPs),
                                                             colors = Colors,animal_labels=Labels_NPHs)

# --- CD40.RBDv-Conv group
Plot_VL_Trachea_CD40RBDvConv <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "CD40.RBDv - Conv" & ObservationType %in% c("gRNA_T") & SubjectID %in% Selected_NHPs),
                                                                     observed_data=subset(ObservedData,Group == "CD40.RBDv - Conv" & MeasureType %in% c("gRNA") & SampleType == "Tracheal fluid" & SubjectID %in% Selected_NHPs),
                                                                     colors = Colors,animal_labels=Labels_NPHs)

# --- CD40.PanCov-Conv group
Plot_VL_Trachea_CD40PanCovConv <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "CD40.PanCoV - Conv" & ObservationType %in% c("gRNA_T") & SubjectID %in% Selected_NHPs),
                                                                       observed_data=subset(ObservedData,Group == "CD40.PanCoV - Conv" & MeasureType %in% c("gRNA") & SampleType == "Tracheal fluid" & SubjectID %in% Selected_NHPs),
                                                                       colors = Colors,animal_labels=Labels_NPHs)

# --- mRNA-Conv group
Plot_VL_Trachea_mRNAConv <- Individual_Plots_gRNA_sgRNA_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "mRNA - Conv" & ObservationType %in% c("gRNA_T") & SubjectID %in% Selected_NHPs),
                                                                 observed_data=subset(ObservedData,Group == "mRNA - Conv" & MeasureType %in% c("gRNA") & SampleType == "Tracheal fluid" & SubjectID %in% Selected_NHPs),
                                                                 colors = Colors,animal_labels=Labels_NPHs)


# --- Merged groups 
Legend <- cowplot::get_legend(Plot_VL_Trachea_Naive)

Figure <- grid.arrange(Plot_VL_Trachea_Naive + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_VL_Trachea_CD40RBDvNaive + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_VL_Trachea_Conv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_VL_Trachea_CD40RBDvConv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_VL_Trachea_CD40PanCovConv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_VL_Trachea_mRNAConv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       layout_matrix = rbind(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3)))

Plot_VL_Trachea_AllGroups <- ggarrange(annotate_figure(Figure, left = textGrob("Viral load (copies/mL)", rot = 90, just="centre", gp = gpar(fontsize = 12,fontface="bold")),
                                                       bottom = textGrob("Time post-exposure (days)",just="centre", gp = gpar(fontsize = 12,fontface="bold")))) + 
  bgcolor("white") + border("white") 
# ggsave(Plot_VL_Trachea_AllGroups,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"IndivFits_VL_Trachea.png",sep="/"),height=12,width=10,dpi=300)



### > Antibody  ---- 
# --- Naive group
Plot_Antibody_Naive <- Individual_Plots_BAb_NAb_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "Naive" & ObservationType %in% c("Obs_BAb","Obs_NAb") & SubjectID %in% Selected_NHPs),
                                                         observed_data=subset(ObservedData,Group == "Naive" & MeasureType %in% c("BAb","NAb")  & SubjectID %in% Selected_NHPs),
                                                         colors = Colors,animal_labels = Labels_NPHs) + guides(fill='none')

# --- CD40.RBDv-Naive Group
Plot_Antibody_CD40RBDvNaive <- Individual_Plots_BAb_NAb_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "CD40.RBDv - Naive" & ObservationType %in% c("Obs_BAb","Obs_NAb")  & SubjectID %in% Selected_NHPs ),
                                                                 observed_data=subset(ObservedData,Group == "CD40.RBDv - Naive" & MeasureType %in% c("BAb","NAb") & SubjectID %in% Selected_NHPs),
                                                                 colors = Colors,animal_labels = Labels_NPHs) + theme(legend.position = "top") + guides(fill='none')

# --- Conv Group
Plot_Antibody_Conv <- Individual_Plots_BAb_NAb_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "Conv" & ObservationType %in% c("Obs_BAb","Obs_NAb") & SubjectID %in% Selected_NHPs),
                                                        observed_data=subset(ObservedData,Group == "Conv" & MeasureType %in% c("BAb","NAb") & SubjectID %in% Selected_NHPs),
                                                        colors = Colors,animal_labels = Labels_NPHs) + guides(fill='none')

# --- CD40.RBDv-Conv group
Plot_Antibody_CD40RBDvConv <- Individual_Plots_BAb_NAb_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "CD40.RBDv - Conv" & ObservationType %in% c("Obs_BAb","Obs_NAb") & SubjectID %in% Selected_NHPs),
                                                                observed_data=subset(ObservedData,Group == "CD40.RBDv - Conv" & MeasureType %in% c("BAb","NAb") & SubjectID %in% Selected_NHPs),
                                                                colors = Colors,animal_labels = Labels_NPHs) + guides(fill='none')

# --- CD40.PanCov-Conv group
Plot_Antibody_CD40PanCovConv <- Individual_Plots_BAb_NAb_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "CD40.PanCoV - Conv" & ObservationType %in% c("Obs_BAb","Obs_NAb") & SubjectID %in% Selected_NHPs),
                                                                  observed_data=subset(ObservedData,Group == "CD40.PanCoV - Conv" & MeasureType %in% c("BAb","NAb") & SubjectID %in% Selected_NHPs),
                                                                  colors = Colors,animal_labels = Labels_NPHs) + guides(fill='none')

# --- mRNA-Conv group
Plot_Antibody_mRNAConv <- Individual_Plots_BAb_NAb_FUNCTION(model_data=subset(Distribution_Individual_dynamics,Group == "mRNA - Conv" & ObservationType %in% c("Obs_BAb","Obs_NAb") & SubjectID %in% Selected_NHPs),
                                                            observed_data=subset(ObservedData,Group == "mRNA - Conv" & MeasureType %in% c("BAb","NAb") & SubjectID %in% Selected_NHPs),
                                                            colors = Colors,animal_labels = Labels_NPHs) + guides(fill='none')

# --- Merged groups 
Legend <- cowplot::get_legend(Plot_Antibody_Naive + theme(legend.direction = "vertical",legend.box = "horizontal",legend.title.position = "left"))

Figure <- grid.arrange(Plot_Antibody_Naive + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_Antibody_CD40RBDvNaive + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_Antibody_Conv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_Antibody_CD40RBDvConv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_Antibody_CD40PanCovConv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       Plot_Antibody_mRNAConv + rremove("ylab") + rremove("xlab") + rremove("legend") + theme(plot.margin = unit(c(0.05,0.05,0.05,0.05), "cm")),
                       layout_matrix = rbind(rep(1,3),rep(2,3),rep(3,3),rep(4,3),rep(5,3),rep(6,3)))


Plot_Antibody_AllGroups <- ggarrange(annotate_figure(Figure,
                                                     left = textGrob("Anti-RBD IgG (AU/mL)", rot = 90, just="centre", gp = gpar(fontsize = 12,fontface="bold")),
                                                     bottom = textGrob("Time post-exposure (days)",just="centre", gp = gpar(fontsize = 12,fontface="bold"))),
                                     legend="top",legend.grob = Legend) + 
  bgcolor("white") + border("white")  
# ggsave(Plot_Antibody_AllGroups,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"IndivFits_Antibody.png",sep="/"),height=13,width=10,dpi=300)
# ------------ #







# ........................----
# -- UNOBSERVED COMPARTMENTS ----
# Paper: Figures S5B ----
# ........................----
## Simulation of population parameters ----
initializeLixoftConnectors(software = "simulx",path = "C:/ProgramData/Lixoft/MonolixSuite2023R1",force = T)
Simulated_Population_Parameters <- simpopmlx(n=Nb_pop,project = Monolix_Project,seed = Seed)
# Removal of inter-individual variability and observation uncertainty
Simulated_Population_Parameters[,which(names(Simulated_Population_Parameters) %in% c(paste("omega",c("delta_N","S0","theta_BAb","eta"),sep="_"),
                                                                                     paste("a",seq(1,6),sep="")))] <- 0



## Study of unobserved compartments ----
### > Simulation of dynamics ----
Time_simulation <- seq(0,30,by=0.25)
Model_Outputs <- list(list(name="TN_T0_Percent",time=Time_simulation),
                      list(name="S",time=Time_simulation))

NHPs_Profile <- data.frame(Group = unique(NHPs$Group),Inoc=median(NHPs$Inoc),Weight=median(NHPs$Weight))
NHPs_Profile <- cbind(id = seq(1,nrow(NHPs_Profile)),NHPs_Profile)

Covariates <- NHPs_Profile[,c("id","Group")]
Regressors <- merge(NHPs_Profile,expand.grid(id=unique(NHPs_Profile$id),time=Time_simulation))
Regressors <- Regressors[order(Regressors$time,Regressors$id),]


Simulated_UnobservedCompartment_Dynamics <- simulx(model = MlxtranModel_Simulation_file,
                                                   parameter = Simulated_Population_Parameters,
                                                   output = Model_Outputs,regressor = Regressors,
                                                   covariate = Covariates,
                                                   settings = list(sharedIds="population"))[c("TN_T0_Percent","S")]
Simulated_UnobservedCompartment_Dynamics <- Reduce(function(x,y) merge(x,y,all=TRUE), Simulated_UnobservedCompartment_Dynamics)
Simulated_UnobservedCompartment_Dynamics <- merge(Covariates,Simulated_UnobservedCompartment_Dynamics)[-1]
names(Simulated_UnobservedCompartment_Dynamics)[which(names(Simulated_UnobservedCompartment_Dynamics) == "rep")] <- "pop"



### > Plot of dynamics ----
# Reshape of the dataset
Simulated_UnobservedCompartment_Dynamics <- melt(data=Simulated_UnobservedCompartment_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")
# Calculation of distributions
Distribution_UnobservedCompartment_Dynamics <- ddply(.data=Simulated_UnobservedCompartment_Dynamics,.variables = .(Group,time,ObservationType),summarize,
                                                     Mean = mean(Value,na.rm=T),
                                                     Median = median(Value,na.rm=T),
                                                     ICMIN = quantile(Value,probs=c(0.025)),
                                                     Q1 = quantile(Value,probs=(0.25)),
                                                     Q3 = quantile(Value,probs=(0.75)),
                                                     ICMAX = quantile(Value,probs=c(0.975)))


Colors <- Plot_colors_Groups()

Plot_S_dynamics <- ggplot(data=subset(Distribution_UnobservedCompartment_Dynamics,ObservationType == "S" & Group %in% c("Naive","Conv"))) + 
  geom_ribbon(aes(x=time,ymin=Q1,ymax=Q3,color=as.factor(Group),fill=as.factor(Group)),linetype="dashed",linewidth=1.0,alpha=0.30) +
  geom_line(aes(x=time,y=Mean,color=as.factor(Group)),linewidth=1.2) + 
  scale_color_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE) + 
  scale_fill_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE)  + 
  
  scale_x_continuous(breaks = seq(0,30,by=5)) + xlab("Time post-exposure (days)") + 
  scale_y_continuous(breaks = seq(0,1000,by=150)) + ylab("Antibody secreting cells (cells)") + 
  coord_cartesian(xlim = c(0,30),ylim=c(-1,800),expand = 0) +
  theme_bw() + 
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=12,face="bold"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        legend.position = "inside",
        legend.position.inside = c(0.15,0.75),
        legend.title = element_blank(),
        legend.key.width = unit(units = "cm",x = 1.5),
        legend.key.height = unit(units = "cm",x = 0.5),
        panel.grid.major = element_line(colour = "gray90"))
# ggsave(Plot_S_dynamics,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"ASCs_dynamics_Naive_Vs_Convalescent.png",sep="/"),height=4,width=6)
# ------------ #







# ........................----
# -- NEUTRALIZATION QUALITY OF ANTIBODIES  ----
# Paper: Figures 6, S5A ----
# ........................----
## Simulation of population parameters ----
initializeLixoftConnectors(software = "simulx",path = "C:/ProgramData/Lixoft/MonolixSuite2023R1",force = T)
Simulated_Population_Parameters <- simpopmlx(n=Nb_pop,project = Monolix_Project,seed = Seed)
# Removal of inter-individual variability and observation uncertainty
Simulated_Population_Parameters[,which(names(Simulated_Population_Parameters) %in% c(paste("omega",c("delta_N","S0","theta_BAb","eta"),sep="_"),
                                                                                     paste("a",seq(1,6),sep="")))] <- 0


## Study of neutralization and antibody dynamics (Population dynamics) ----
### > Simulation of dynamics ----
Time_simulation <- seq(0,30,by=0.25)
Model_Outputs <- list(list(name="epsilon",time=Time_simulation),
                      list(name="NAb",time=Time_simulation))

NHPs_Profile <- data.frame(Group = unique(NHPs$Group),Inoc=median(NHPs$Inoc),Weight=median(NHPs$Weight))
NHPs_Profile <- cbind(id = seq(1,nrow(NHPs_Profile)),NHPs_Profile)

Covariates <- NHPs_Profile[,c("id","Group")]
Regressors <- merge(NHPs_Profile,expand.grid(id=unique(NHPs_Profile$id),time=Time_simulation))
Regressors <- Regressors[order(Regressors$time,Regressors$id),]


Simulated_Neutralization_Dynamics <- simulx(model = MlxtranModel_Simulation_file,
                                            parameter = Simulated_Population_Parameters,
                                            output = Model_Outputs,regressor = Regressors,
                                            covariate = Covariates,
                                            settings = list(sharedIds="population"))[c("NAb","epsilon")]
Simulated_Neutralization_Dynamics <- Reduce(function(x,y) merge(x,y,all=TRUE), Simulated_Neutralization_Dynamics)
Simulated_Neutralization_Dynamics <- merge(Covariates,Simulated_Neutralization_Dynamics)[-1]
names(Simulated_Neutralization_Dynamics)[which(names(Simulated_Neutralization_Dynamics) == "rep")] <- "pop"


### > Plot of dynamics (Paper, Figure 6A)  ----
Simulated_Neutralization_Dynamics$epsilon <- Simulated_Neutralization_Dynamics$epsilon*100
# Reshape of the dataset
Simulated_Neutralization_Dynamics <- melt(data=Simulated_Neutralization_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")
# Calculation of distributions
Distribution_Neutralization_Dynamics <- ddply(.data=Simulated_Neutralization_Dynamics,.variables = .(Group,time,ObservationType),summarize,
                                              Mean = mean(Value,na.rm=T),
                                              Median = median(Value,na.rm=T),
                                              ICMIN = quantile(Value,probs=c(0.025)),
                                              Q1 = quantile(Value,probs=(0.25)),
                                              Q3 = quantile(Value,probs=(0.75)),
                                              ICMAX = quantile(Value,probs=c(0.975)))

Colors <- Plot_colors_Groups()

Plot_NAb_dynamics <- ggplot(data=subset(Distribution_Neutralization_Dynamics,ObservationType == "NAb")) + 
  geom_hline(yintercept = 0.78,linetype="solid",color="darkred",linewidth=1.0,alpha=0.8) +
  geom_ribbon(aes(x=time,ymin=ICMIN,ymax=ICMAX,color=as.factor(Group),fill=as.factor(Group)),linetype="dashed",linewidth=1.0,alpha=0.30) + 
  geom_line(aes(x=time,y=Mean,color=as.factor(Group)),linewidth=1.2) + 
  scale_color_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE) + 
  scale_fill_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE)  + 
  
  scale_x_continuous(breaks = seq(0,30,by=5)) + xlab("Time post-exposure (days)") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = exp_bold) + ylab("ACE2/RBD binding inhibition (AU/mL)") + 
  coord_cartesian(xlim = c(0,30),ylim=c(0.5e-1,1E3),expand = 0) +
  
  theme_bw() + 
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=12,face="bold"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.key.width = unit(units = "cm",x = 1.5),
        legend.key.height = unit(units = "cm",x = 0.5),
        panel.grid.major = element_line(colour = "gray90"))

data.tb <- tibble(x = 29.8, y = 1,plot = list(Plot_NAb_dynamics + ylab("ACE2/RBD bind. inhibit° \n (AU/mL)") +  
                                                ggtitle("ACE2/RBD bind. inhibit° (AU/mL)") + 
                                                theme(legend.position="none",
                                                      panel.grid = element_blank(),
                                                      panel.grid.major = element_blank(),
                                                      plot.title = element_text(color="black",face="bold",size=10),
                                                      axis.title = element_blank(),
                                                      plot.background = element_rect(linewidth = 1.2,fill = alpha("white",1),color="darkblue"))
))

Plot_Epsilon_dynamics <- ggplot(data=subset(Distribution_Neutralization_Dynamics,ObservationType == "epsilon")) + 
  geom_ribbon(aes(x=time,ymin=ICMIN,ymax=ICMAX,color=as.factor(Group),fill=as.factor(Group)),linetype="dashed",linewidth=1.0,alpha=0.30) + 
  geom_line(aes(x=time,y=Mean,color=as.factor(Group)),linewidth=1.2) + 
  
  scale_color_manual(name=NULL,breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE) + 
  scale_fill_manual(name=NULL,breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE)  + 
  
  scale_x_continuous(breaks = seq(0,30,by=5)) + xlab("Time post-exposure (days)") + 
  scale_y_continuous(breaks = seq(0,100,by=20)) + ylab("Neutralization Efficacy (%)") + 
  coord_cartesian(xlim = c(0,30),ylim=c(0,100),expand = 0) +
  
  # Inclusion of the plot of NAb dynamics
  geom_plot(data=data.tb,aes(x,y,label=plot,vp.height = 0.50,vp.width = 0.65)) + 
  guides(color=guide_legend(nrow=2,override.aes = list(linewidth=0.80))) + 
  
  theme_bw() + 
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=11,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=13,face="bold"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        legend.text = element_text(color="black",face="bold",size = 10),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.box = "horizontal",
        legend.title = element_text(face = "bold"),
        legend.key.width = unit(units = "cm",x = 1.5),
        legend.key.height = unit(units = "cm",x = 0.5),
        panel.grid.major = element_line(colour = "gray90"))
# ggsave(Plot_Epsilon_dynamics,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"NeutralizationEfficacy_dynamics.png",sep="/"),height=6,width=10,dpi=300)




## Study of neutralization efficacy vs antibody levels ----
### > Simulation of dynamics ----
Tested_NAb <- 10^seq(-2,3,by=0.1)
Simulated_Population_Parameters_Group <- do.call("rbind",lapply(1:Nb_pop,function(l) rbind(cbind(pop=l,URT_Comp="Naso",Extraction_Parameters_Compartment_FUNCTION(parameters=Simulated_Population_Parameters[l,])$Nasopharynx),
                                                                                           cbind(pop=l,URT_Comp="Trachea",Extraction_Parameters_Compartment_FUNCTION(parameters=Simulated_Population_Parameters[l,])$Trachea))))

Epsilon_Vs_NAb <- NULL
for(i in 1:nrow(Simulated_Population_Parameters_Group)){
  # i <- 1
  if(i%%100 == 0){print(i)}
  
  params <- Simulated_Population_Parameters_Group[i,]
  simulated_res <- data.frame(params[,c("pop","URT_Comp","Group","eta")],NAb=Tested_NAb)
  simulated_res$epsilon <- sapply(1:nrow(simulated_res),function(k) Epsilon_FUNCTION(eta=simulated_res$eta[k],NAb = simulated_res$NAb[k]))
  
  Epsilon_Vs_NAb <- rbind(Epsilon_Vs_NAb,simulated_res)
}


### > Plot of dynamics (Paper, Figure 6B) ----
# library(gghighlight)
Summary_Epsilon_Vs_NAb <- ddply(.data=Epsilon_Vs_NAb,.variables = .(URT_Comp,Group,NAb),summarise,
                                Mean=mean(100*epsilon,na.rm=TRUE),
                                Median=median(100*epsilon,na.rm=TRUE),
                                ICMIN=quantile(100*epsilon,probs=c(0.025)),
                                ICMAX=quantile(100*epsilon,probs=c(0.975)))

Summary_Epsilon_Vs_NAb$Group <- as.character(Summary_Epsilon_Vs_NAb$Group)
Summary_Epsilon_Vs_NAb$Group2 <- with(Summary_Epsilon_Vs_NAb,ifelse(Group %in% c("Naive","CD40.RBDv - Naive","Conv"),"No - ConvXVacc",Group))
Summary_Epsilon_Vs_NAb$Group2 <- factor(Summary_Epsilon_Vs_NAb$Group2,levels= c("No - ConvXVacc","CD40.RBDv - Conv","mRNA - Conv","CD40.PanCoV - Conv")) 

Labels_Group <- c("Naive | CD40.RBDv - Naive | Conv","CD40.RBDv - Conv","mRNA - Conv","CD40.PanCoV - Conv")
names(Labels_Group) <- c("No - ConvXVacc","CD40.RBDv - Conv","mRNA - Conv","CD40.PanCoV - Conv")

Strip_settings <- strip_themed(background_x = elem_list_rect(fill=alpha(c(rgb(125,148,136,maxColorValue = 255),Colors[4:6]),alpha = 0.25),
                                                             color=alpha(c(rgb(125,148,136,maxColorValue = 255),Colors[4:6]),alpha = 1)),
                               text_x = list(element_text(size = 11,face="bold",margin = margin(b=1,t=1)),element_text(size = 11,face="bold",margin = margin(b=1,t=1)),
                                             element_text(size = 11,face="bold",margin = margin(b=1,t=1)),element_text(size = 11,face="bold",margin = margin(b=1,t=1))),
                               by_layer_x = FALSE)


Plot_Epsilon_Vs_NAb <- ggplot(data=subset(Summary_Epsilon_Vs_NAb,URT_Comp == "Naso")) +
  geom_line(aes(x=NAb,y=Mean,color=as.factor(Group)),linetype="solid",linewidth=1.2) + 
  geom_ribbon(aes(x=NAb,ymin=ICMIN,ymax=ICMAX,fill=as.factor(Group)),alpha=0.2,linewidth=1.0,linetype="dashed") +
  geom_ribbon(aes(x=NAb,ymin=ICMIN,ymax=ICMAX,color=as.factor(Group)),alpha=1,linetype="dashed",fill=NA,linewidth=1.0) +
  
  facet_wrap2(~Group2,labeller = labeller(Group2 = Labels_Group),strip=Strip_settings) + 
  
  scale_fill_manual(name="Groups",breaks = names(Colors),values = Colors) + 
  scale_color_manual(name="Groups",breaks = names(Colors),values = Colors) + 
  
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = exp_bold) + 
  xlab("ACE2/RBD binding inhibition (AU/mL)") + 
  
  ylab("Neutralization Efficacy (%)") + 
  coord_cartesian(ylim=c(-1,101),xlim=c(1e-2,1.1e3),expand = 0) +
  theme_bw()  + 
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text.x =  element_text(color="black",size=11,face="bold",hjust = 0.8),
        axis.text.y =  element_text(color="black",size=11,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=13,face="bold"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        legend.position = "none",
        panel.spacing.x = unit(0.75,units="cm"),
        panel.grid.major = element_line(colour = "gray90")#,
        # strip.text = element_text(size = 11,face="bold")
  ) # + 
# gghighlight(unhighlighted_params = list(fill=NA,linewidth=0.70))
# ggsave(Plot_Epsilon_Vs_NAb,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"Neutralization_Vs_Antibody.png",sep="/"),height=6,width=10,dpi=300)



### > Merged plots (Paper, Figure 6) ----
Merged_plot <- plot_grid(Plot_Epsilon_dynamics,Plot_Epsilon_Vs_NAb,
                         labels = c("A","B"),nrow=2,rel_heights = c(2,2))
# ggsave(Merged_plot,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"Study_NeutralizationEfficacy.png",sep="/"),height=12,width=10,dpi=300)




## Study of Effective viral infectivity over time ----
### > Simulation of dynamics ----
Time_simulation <- seq(0,30,by=0.25)
Model_Outputs <- list(list(name="Eff_beta_N",time=Time_simulation))

NHPs_Profile <- data.frame(Group = unique(NHPs$Group),Inoc=median(NHPs$Inoc),Weight=median(NHPs$Weight))
NHPs_Profile <- cbind(id = seq(1,nrow(NHPs_Profile)),NHPs_Profile)

Covariates <- NHPs_Profile[,c("id","Group")]
Regressors <- merge(NHPs_Profile,expand.grid(id=unique(NHPs_Profile$id),time=Time_simulation))
Regressors <- Regressors[order(Regressors$time,Regressors$id),]


Simulated_Infectivity_Dynamics <- simulx(model = MlxtranModel_Simulation_file,
                                         parameter = Simulated_Population_Parameters,
                                         output = Model_Outputs,regressor = Regressors,
                                         covariate = Covariates,
                                         settings = list(sharedIds="population"))[["Eff_beta_N"]]

Simulated_Infectivity_Dynamics <- merge(Covariates,Simulated_Infectivity_Dynamics)[-1]
names(Simulated_Infectivity_Dynamics)[which(names(Simulated_Infectivity_Dynamics) == "rep")] <- "pop"


### > Plot of dynamics (Paper, Figure S5A) ----
# Calculation of distributions
Distribution_Infectivity_Dynamics <- ddply(.data=Simulated_Infectivity_Dynamics,.variables = .(Group,time),summarize,
                                           Mean = mean(Eff_beta_N,na.rm=T),
                                           Median = median(Eff_beta_N,na.rm=T),
                                           ICMIN = quantile(Eff_beta_N,probs=c(0.05)),
                                           Q1 = quantile(Eff_beta_N,probs=(0.25)),
                                           Q3 = quantile(Eff_beta_N,probs=(0.75)),
                                           ICMAX = quantile(Eff_beta_N,probs=c(0.95)))

Colors <- Plot_colors_Groups()
Distribution_Infectivity_Dynamics$Group <- factor(Distribution_Infectivity_Dynamics$Group,levels = names(Colors))

Plot_Infectivity_dynamics <- ggplot(data=Distribution_Infectivity_Dynamics) + 
  geom_ribbon(aes(x=time,ymin=ICMIN,ymax=ICMAX,color=as.factor(Group),fill=as.factor(Group)),linetype="dashed",linewidth=1.0,alpha=0.30) + 
  geom_line(aes(x=time,y=Mean,color=as.factor(Group)),linewidth=1.2) + 
  
  scale_color_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE) + 
  scale_fill_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE)  + 
  
  scale_x_continuous(breaks = seq(0,30,by=5)) + xlab("Time post-exposure (days)") + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = exp_bold) + 
  ylab("Effective Viral infectivity (ml/(copies.day))") + 
  coord_cartesian(xlim = c(0,30), ylim=c(1.5E-10,0.25E-3),expand = 0) +
  
  guides(color=guide_legend(nrow=1)) + 
  theme_bw() + 
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=12,face="bold"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        legend.position = "top",
        legend.direction = "horizontal",
        legend.title = element_text(face = "bold"),
        legend.key.width = unit(units = "cm",x = 1.5),
        legend.key.height = unit(units = "cm",x = 0.5),
        panel.grid.major = element_line(colour = "gray90"))
# ggsave(Plot_Infectivity_dynamics,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"Effective_Viral_Infectivity.png",sep="/"),height=5,width=11)




### > Merged plot with ASCs dynamics (Paper Figure S5) ----
# Merged plot generated for the modeling paper
# To obtain this merged plot, the section "Unobserved Compartments must be launched before to get the plot "Plot_S_dynamics"

# We use the global legend of the plot "Plot_Infectivity_dynamics"
Legend <- ggpubr::get_legend(Plot_Infectivity_dynamics)
Merged_Plot_Infectivity_ASCs_without_legend <- plot_grid(Plot_Infectivity_dynamics + theme(legend.position = "none") + coord_cartesian(xlim = c(0,25), ylim=c(1.5E-10,0.25E-3),expand = 0),
                                                         Plot_S_dynamics + theme(legend.position = "none") + coord_cartesian(xlim = c(0,25),ylim=c(-1,800),expand = 0),
                                                         labels = c("A","B"),ncol=2)
Merged_Plot_Infectivity_ASCs <- plot_grid(Legend,Merged_Plot_Infectivity_ASCs_without_legend,
                                          ncol=1,rel_heights = c(0.05,1))
# ggsave(Merged_Plot_Infectivity_ASCs,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"Effective_Viral_Infectivity_and_ASCs_Dynamics.png",sep="/"),height=5,width=10,bg = "white")






## Study of Basic reproduction number vs antibodies levels ----
### > Simulation of dynamics ----
Tested_NAb <- 10^seq(-2,3,by=0.1)
Simulated_Population_Parameters_Group <- do.call("rbind",lapply(1:Nb_pop,function(l) rbind(cbind(pop=l,URT_Comp="Naso",Extraction_Parameters_Compartment_FUNCTION(parameters=Simulated_Population_Parameters[l,])$Nasopharynx),
                                                                                           cbind(pop=l,URT_Comp="Trachea",Extraction_Parameters_Compartment_FUNCTION(parameters=Simulated_Population_Parameters[l,])$Trachea))))
NHPs_Profile <- data.frame(Group = unique(NHPs$Group),Inoc=median(NHPs$Inoc),Weight=median(NHPs$Weight))
NHPs_Profile <- cbind(id = seq(1,nrow(NHPs_Profile)),NHPs_Profile)

R0_Vs_NAb <- NULL
for(i in 1:Nb_pop){
  # i <- 1
  if(i%%100 == 0){print(i)}
  
  params <- subset(Simulated_Population_Parameters_Group,pop == i)
  
  for(group in unique(params$Group)){
    # group <- unique(params$Group)[1]
    
    params_group <- subset(params,Group == group)
    combined_params_group <- unique(params_group[,names(params_group) %notin% c("URT_Comp","P")])
    combined_params_group$P_N <- subset(params_group,URT_Comp == "Naso")$P
    combined_params_group$P_T <- subset(params_group,URT_Comp == "Trachea")$P
    combined_params_group$Weight <- subset(NHPs_Profile,Group == combined_params_group$Group)$Weight
    combined_params_group$WT <- with(combined_params_group,ifelse(Weight<=4.5,2,3))
    combined_params_group$WN <- with(combined_params_group,ifelse(Weight<=4.5,4,5.5))
    combined_params_group$T0_T <- 2.25*10^4/combined_params_group$WT
    combined_params_group$T0_N <- 1.25*10^5/combined_params_group$WN   
    
    simulated_res <- data.frame(combined_params_group,NAb_0=Tested_NAb)
    simulated_res$R0 <- sapply(1:nrow(simulated_res),function(k) Basic_Reproduction_Number_FUNCTION(parameters = simulated_res[k,]))
    
    R0_Vs_NAb <- rbind(R0_Vs_NAb,simulated_res[,c("pop","Group","NAb_0","R0")])
  }
}

### > Plot of dynamics ----
library(gghighlight)

Distribution_R0_Vs_NAb <- ddply(.data=R0_Vs_NAb,.variables = .(Group,NAb_0),summarise,
                                Mean=mean(R0,na.rm=TRUE),
                                Median=median(R0,na.rm=TRUE),
                                ICMIN=quantile(R0,probs=c(0.025)),
                                ICMAX=quantile(R0,probs=c(0.975)))

Distribution_R0_Vs_NAb$Group <- factor(Distribution_R0_Vs_NAb$Group,levels = names(Colors))
Distribution_R0_Vs_NAb$Group2 <- factor(Distribution_R0_Vs_NAb$Group,levels = names(Colors))

Strip_settings <- strip_themed(background_x = elem_list_rect(fill=alpha(Colors,alpha = 0.25),
                                                             color=alpha(Colors,alpha = 1)),
                               by_layer_x = FALSE)

Plot_R0_Vs_NAb <- ggplot(data=Distribution_R0_Vs_NAb) +
  geom_hline(yintercept = 1,color="darkred",linewidth=1.0,linetype="solid") +
  
  geom_line(aes(x=NAb_0,y=Mean,color=as.factor(Group)),linetype="solid",linewidth=1.2) + 
  geom_ribbon(aes(x=NAb_0,ymin=ICMIN,ymax=ICMAX,fill=as.factor(Group)),alpha=0.2,linewidth=1.0) +
  geom_ribbon(aes(x=NAb_0,ymin=ICMIN,ymax=ICMAX,color=as.factor(Group)),alpha=1,linetype="dashed",fill=NA,linewidth=1.0) +
  
  facet_wrap2(~Group2,strip=Strip_settings) + 
  
  scale_fill_manual(name="Groups",breaks = names(Colors),values = Colors) + 
  scale_color_manual(name="Groups",breaks = names(Colors),values = Colors) + 
  
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = exp_bold) + 
  xlab("ACE2/RBD binding inhibition (AU/mL)") + 
  
  ylab("Basic Reproduction Number (R0)") + 
  coord_cartesian(ylim=c(0,9),xlim=range(Tested_NAb),expand = 0) +
  theme_bw()  + 
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=12,face="bold"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        legend.position = "none",
        panel.spacing.x = unit(0.75,units="cm"),
        panel.grid.major = element_line(colour = "gray90"),
        strip.text = element_text(size = 9,face="bold")) + 
  gghighlight(unhighlighted_params = list(fill=NA,linewidth=0.70))

# ggsave(Plot_R0_Vs_NAb,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"BasicReproduction_Vs_Antibody.png",sep="/"),height=6,width=12)
# ------------ #





# ........................----
# -- PROTECTIVE THRESHOLD AND BASIC REPRODUCTION NUMBER ----
# Paper: Figure 7 ----
# ........................----
initializeLixoftConnectors(software = "simulx",path = "C:/ProgramData/Lixoft/MonolixSuite2023R1",force = T)
Simulated_Population_Parameters <- simpopmlx(n=Nb_pop,project = Monolix_Project,seed = Seed)
# Removal of inter-individual variability and observation uncertainty
Simulated_Population_Parameters[,which(names(Simulated_Population_Parameters) %in% c(paste("omega",c("delta_N","S0","theta_BAb","eta"),sep="_"),
                                                                                     paste("a",seq(1,6),sep="")))] <- 0


## Study of protective threshold and R0 ----
### > Simulation of results ----
Time_simulation <- c(0)
Model_Outputs <- list(list(name="R0",time=Time_simulation),
                      list(name="NAb_Thresh",time=Time_simulation))

NHPs_Profile <- data.frame(Group = unique(NHPs$Group),Inoc=median(NHPs$Inoc),Weight=median(NHPs$Weight))
NHPs_Profile <- cbind(id = seq(1,nrow(NHPs_Profile)),NHPs_Profile)

Covariates <- NHPs_Profile[,c("id","Group")]
Regressors <- merge(NHPs_Profile,expand.grid(id=unique(NHPs_Profile$id),time=Time_simulation))
Regressors <- Regressors[order(Regressors$time,Regressors$id),]

Simulated_Protective_Threshold <- simulx(model = MlxtranModel_Simulation_file,
                                         parameter = Simulated_Population_Parameters,
                                         output = Model_Outputs,regressor = Regressors,
                                         covariate = Covariates,
                                         settings = list(sharedIds="population"))[c("NAb_Thresh","R0")]

Simulated_Protective_Threshold <- Reduce(function(x,y) merge(x,y,all=TRUE), Simulated_Protective_Threshold)
Simulated_Protective_Threshold <- merge(Covariates,Simulated_Protective_Threshold)
Simulated_Protective_Threshold <- Simulated_Protective_Threshold[,!names(Simulated_Protective_Threshold) %in% c("id","time")]
names(Simulated_Protective_Threshold)[which(names(Simulated_Protective_Threshold) == "rep")] <- "pop"

### > Plot of results (Figures 7A, 7B) ----
library(ggrepel)
# Reshape of the dataset
Simulated_Protective_Threshold <- melt(data=Simulated_Protective_Threshold,id.vars = c("pop","Group"),variable.name = "ObservationType",value.name = "Value")
# Calculation of distributions
Distribution_Protective_Thresold <- ddply(.data=Simulated_Protective_Threshold,.variables = .(Group,ObservationType),summarize,
                                          Mean = mean(Value,na.rm=T),
                                          Median = median(Value,na.rm=T),
                                          ICMIN = quantile(Value,probs=c(0.025)),
                                          Q1 = quantile(Value,probs=(0.25)),
                                          Q3 = quantile(Value,probs=(0.75)),
                                          ICMAX = quantile(Value,probs=c(0.975)))

Distribution_Protective_Thresold$Label <- with(Distribution_Protective_Thresold,
                                               paste(format(Mean,digits=1,scientific=TRUE), " [",format(ICMIN,digits=2,scientific=TRUE)," ; ",format(ICMAX,digits=2,scientific=TRUE),"]",sep=""))


Colors <- Plot_colors_Groups()
Simulated_Protective_Threshold$Group <- factor(Simulated_Protective_Threshold$Group,levels = names(Colors))
Distribution_Protective_Thresold$Group <- factor(Distribution_Protective_Thresold$Group,levels = names(Colors))
ObservedData$Group <- factor(ObservedData$Group,levels = names(Colors))

Strip_settings <- strip_themed(background_x = elem_list_rect(fill=alpha(Colors,alpha = 0.25),color=alpha(Colors,alpha = 1)),by_layer_x = FALSE)

Plot_Protective_Threshold <- ggplot(data=subset(Distribution_Protective_Thresold,ObservationType == "NAb_Thresh")) + 
  geom_boxplot(aes(x=as.factor(Group),ymin=ICMIN,lower=Q1,middle=Median,upper=Q3,ymax=ICMAX,fill=as.factor(Group)),stat="identity",color="black",linewidth=1.0,alpha=0.5,outlier.shape = NA) + 
  geom_point(data=subset(Distribution_Protective_Thresold,ObservationType == "NAb_Thresh"), aes(x=Group,y=ICMAX), pch=16, cex=2,color="red") +
  
  geom_label_repel(data=subset(Distribution_Protective_Thresold,ObservationType == "NAb_Thresh"),aes(x=Group,y=ICMAX,label=round(ICMAX,digits=1),color=as.factor(Group),fill=as.factor(Group)),
                   force = 10,nudge_x = -0.25,nudge_y = 0.30,label.size = 1,size=4,alpha=0.20) + 
  geom_label_repel(data=subset(Distribution_Protective_Thresold,ObservationType == "NAb_Thresh"),aes(x=Group,y=ICMAX,label=round(ICMAX,digits=1)),
                   force = 10,nudge_x = -0.25,nudge_y = 0.30,label.size = 1,size=4,alpha=1,fill=NA,color="black") + 
  
  geom_point(data=subset(ObservedData,obsid == 6  & DayPostExpo %in% c(0,3,9)),
             aes(x=Group,y=CensoredValue,fill=as.factor(Group),group=DayPostExpo,shape=as.factor(DayPostExpo)),
             cex=2,stroke=1.5,alpha=0.5,position = position_dodge(width = 1.0)) +
  geom_point(data=subset(ObservedData,obsid == 6  & DayPostExpo %in% c(0,3,9)),
             aes(x=Group,y=CensoredValue,group=DayPostExpo,shape=as.factor(DayPostExpo)),
             color="black",fill=NA,cex=2,stroke=1.5,alpha=1,position = position_dodge(width = 1.0)) +
  
  facet_wrap2(.~Group,nrow=1,scales = "free_x",strip=Strip_settings) + 
  
  scale_fill_manual(name="Groups",breaks = names(Colors),values = Colors) + 
  scale_color_manual(name="Groups",breaks = names(Colors),values = Colors) + 
  scale_shape_manual(name="Observations",breaks=c(0,3,9),labels=c("0 dpe","3 dpe","9 dpe"),values=c(21,24,22)) + 
  guides(color="none",fill="none") + 
  
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = exp_bold) + 
  ylab("Protective threshold (AU/mL)") + 
  
  coord_cartesian(xlim=c(0.25,1.75),ylim = c(2e-3,1e3),expand = 0) +
  theme_bw() +
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=12,face="bold"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.text = element_text(color="black",size=11),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        legend.key.width = unit(units = "cm",x = 0.5),
        legend.position = "top",
        panel.spacing.x = unit(0,units="cm"),
        panel.grid.major = element_line(colour = "gray90"),
        strip.text = element_text(size = 9,face="bold")) 
# ggsave(Plot_Protective_Threshold,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"Protective_Threshold.png",sep="/"),height=4,width=12,dpi=300)


Plot_BasicReproductionNumber <- ggplot(data=subset(Distribution_Protective_Thresold,ObservationType == "R0")) + 
  geom_hline(yintercept = 1,color="red",linewidth=1.0,linetype="dashed") + 
  geom_boxplot(aes(x=as.factor(Group),ymin=ICMIN,lower=Q1,middle=Median,upper=Q3,ymax=ICMAX,fill=as.factor(Group)),stat="identity",color="black",linewidth=1.0,alpha=0.5,outlier.shape = NA) + 
  
  facet_wrap2(.~Group,nrow=1,scales = "free_x",strip=Strip_settings) + 
  
  scale_fill_manual(name="Groups",breaks = names(Colors),values = Colors) + 
  scale_color_manual(name="Groups",breaks = names(Colors),values = Colors) + 
  scale_shape_manual(name="Observations",breaks=c(0,3,9),labels=c("0 dpe","3 dpe","9 dpe"),values=c(21,24,22)) + 
  guides(color="none",fill="none") + 
  
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = exp_bold) + 
  ylab("Basic reproduction number") + 
  
  coord_cartesian(xlim=c(0.25,1.75),ylim = c(0.4e-4,10),expand = 0) +
  theme_bw() +
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=12,face="bold"),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.title = element_text(face = "bold"),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        legend.key.width = unit(units = "cm",x = 0.5),
        legend.position = "top",
        panel.spacing.x = unit(0,units="cm"),
        panel.grid.major = element_line(colour = "gray90"),
        strip.text = element_text(size = 9,face="bold")) 
# ggsave(Plot_BasicReproductionNumber,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"Basic_Reproduction_Number.png",sep="/"),height=3.4,width=12,dpi=300)


### > Conversion of Nuetralization threshold into WHO BAU/mL standard ----
AU_to_BAU_Factor <- 0.0272
Extracted_Simulated_Protective_Threshold <- subset(Simulated_Protective_Threshold,ObservationType == "NAb_Thresh")
Extracted_Simulated_alphaNAb <- Simulated_Population_Parameters[,c("pop","alpha_NAb_pop")]
Extracted_Simulated_Protective_Threshold <- merge(Extracted_Simulated_Protective_Threshold,Extracted_Simulated_alphaNAb)
Extracted_Simulated_Protective_Threshold$BAb_Thresh <- with(Extracted_Simulated_Protective_Threshold,Value/alpha_NAb_pop)
Extracted_Simulated_Protective_Threshold$BAU_Thresh <- Extracted_Simulated_Protective_Threshold$BAb_Thresh*AU_to_BAU_Factor

Distribution_BAU_Threshold <-  ddply(.data=Extracted_Simulated_Protective_Threshold,.variables = .(Group),summarize,
                                     Mean = mean(BAU_Thresh,na.rm=T),
                                     Median = median(BAU_Thresh,na.rm=T),
                                     ICMIN = quantile(BAU_Thresh,probs=c(0.025)),
                                     Q1 = quantile(BAU_Thresh,probs=(0.25)),
                                     Q3 = quantile(BAU_Thresh,probs=(0.75)),
                                     ICMAX = quantile(BAU_Thresh,probs=c(0.975)))

Distribution_BAU_Threshold$Label.Mean <- with(Distribution_BAU_Threshold,
                                              paste(format(Mean,digits=2,scientific=FALSE), " [",format(ICMIN,digits=2,scientific=FALSE)," ; ",format(ICMAX,digits=2,scientific=FALSE),"]",sep=""))
Distribution_BAU_Threshold$Label.Median <- with(Distribution_BAU_Threshold,
                                                paste(format(Median,digits=2,scientific=FALSE), " [",format(Q1,digits=2,scientific=FALSE)," ; ",format(Q3,digits=2,scientific=FALSE),"]",sep=""))


## Study of reproduction number dependent on time ----
### > Simulation of results ----
Time_simulation <- seq(0,30,by=0.25)
Model_Outputs <- list(list(name="R0t",time=Time_simulation),
                      list(name="Rt",time=Time_simulation))

NHPs_Profile <- data.frame(Group = unique(NHPs$Group),Inoc=median(NHPs$Inoc),Weight=median(NHPs$Weight))
NHPs_Profile <- cbind(id = seq(1,nrow(NHPs_Profile)),NHPs_Profile)

Covariates <- NHPs_Profile[,c("id","Group")]
Regressors <- merge(NHPs_Profile,expand.grid(id=unique(NHPs_Profile$id),time=Time_simulation))
Regressors <- Regressors[order(Regressors$time,Regressors$id),]

Simulated_ReproductionNumber_Dynamics <- simulx(model = MlxtranModel_Simulation_file,
                                                parameter = Simulated_Population_Parameters,
                                                output = Model_Outputs,regressor = Regressors,
                                                covariate = Covariates,
                                                settings = list(sharedIds="population"))[c("R0t","Rt")]

Simulated_ReproductionNumber_Dynamics <- Reduce(function(x,y) merge(x,y,all=TRUE), Simulated_ReproductionNumber_Dynamics)
Simulated_ReproductionNumber_Dynamics <- merge(Covariates,Simulated_ReproductionNumber_Dynamics)[-1]
names(Simulated_ReproductionNumber_Dynamics)[which(names(Simulated_ReproductionNumber_Dynamics) == "rep")] <- "pop"


### > Plot of results (Figure 7C) ----
# Reshape of the dataset
Simulated_ReproductionNumber_Dynamics <- melt(data=Simulated_ReproductionNumber_Dynamics,id.vars = c("pop","Group","time"),variable.name = "ObservationType",value.name = "Value")
# Calculation of distributions
Distribution_ReproductionNumber <- ddply(.data=Simulated_ReproductionNumber_Dynamics,.variables = .(Group,time,ObservationType),summarize,
                                         Mean = mean(Value,na.rm=T),
                                         Median = median(Value,na.rm=T),
                                         ICMIN = quantile(Value,probs=c(0.025)),
                                         Q1 = quantile(Value,probs=(0.25)),
                                         Q3 = quantile(Value,probs=(0.75)),
                                         ICMAX = quantile(Value,probs=c(0.975)))

Colors <- Plot_colors_Groups()
Distribution_ReproductionNumber$Group <- factor(Distribution_ReproductionNumber$Group,levels = names(Colors))

Strip_settings <- strip_themed(background_x = elem_list_rect(fill=alpha(Colors,alpha = 0.25),color=alpha(Colors,alpha = 1)),by_layer_x = FALSE)

Plot_Reproduction_Number <- ggplot(data=Distribution_ReproductionNumber) + 
  geom_hline(yintercept = 1,color="red",linetype="dashed",linewidth=1.0) + 
  # Confidence Interval Rt
  geom_ribbon(data=subset(Distribution_ReproductionNumber,ObservationType=="Rt"),
              aes(x=time,ymin=ICMIN,ymax=ICMAX,color=as.factor(Group)),linetype="longdash",fill="white",linewidth=0.7,alpha=0.3) +
  geom_ribbon(data=subset(Distribution_ReproductionNumber,ObservationType=="Rt"),
              aes(x=time,ymin=ICMIN,ymax=ICMAX,color=as.factor(Group)),linetype="longdash",fill="white",linewidth=0.7,alpha=1) +
  # Confidence Interval R0t
  geom_ribbon(data=subset(Distribution_ReproductionNumber,ObservationType=="R0t"),
              aes(x=time,ymin=ICMIN,ymax=ICMAX,fill=as.factor(Group),color=as.factor(Group)),linetype="solid",linewidth=0.7,alpha=0.3) +
  geom_ribbon(data=subset(Distribution_ReproductionNumber,ObservationType=="R0t"),
              aes(x=time,ymin=ICMIN,ymax=ICMAX,color=as.factor(Group)),linetype="solid",fill=NA,linewidth=0.7,alpha=1) +
  
  geom_line(aes(x=time,y=Mean,color=as.factor(Group),linetype=as.factor(ObservationType)),linewidth=1.2) + 
  
  facet_wrap2(.~Group,nrow=1,scales = "free_x",strip=Strip_settings) +
  
  scale_color_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE) +
  scale_fill_manual(name="Group",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE)  +
  scale_linetype_manual(name="Reproduction Number",breaks = c("R0t","Rt"),values=c("solid","longdash"),labels=c("Without Target cell depletion","With Target cell depletion")) + 
  guides(fill = 'none', color='none',linetype=guide_legend(nrow=1)) + 
  
  scale_x_continuous(breaks = seq(0,30,by=5)) + xlab("Time post-exposure (days)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = exp_bold) + ylab("Reproduction number") +
  coord_cartesian(xlim = c(0,29.5),ylim=c(0.5e-4,10),expand=0) +
  
  theme_bw() +
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=12,face="bold"),
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
        strip.text = element_text(size = 9,face="bold")) 
# ggsave(Plot_Reproduction_Number,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"Reproduction_Number.png",sep="/"),height=3.5,width=12,dpi=300)


## Merged plots (Figure 7) ----
Merged_Plot_Threshod_and_R0 <- plot_grid(Plot_Protective_Threshold,Plot_BasicReproductionNumber,Plot_Reproduction_Number,
                                         labels=c("A","B","C"),label_fontface = "bold",ncol=1,rel_heights = c(1,0.80,1.0))
# ggsave(Merged_Plot_Threshod_and_R0,filename=paste(MonolixModel_Folder,MonolixModel_name,Figure_Folder,"MergedPlots_ProtectiveThreshold.png",sep="/"),height=10,width=10,dpi = 300)
# ------------ #