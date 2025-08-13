# --------------- # 
# DESCRIPTION: Descriptive analysis of antibody dynamics used in the modeling work
#
# Author: Marie Alexandre 
# R version: 4.2.1

# Open the document outline (Ctrl + shift + O) to see the structure of the document
# --------------- # 




rm(list=ls())
'%notin%' <- Negate('%in%')




# --- LIBRAIRIES ----
library(dplyr)  # library for the function 'case_when' 
library(plyr)   # library to summarize data 
library(reshape2)
library(ggplot2) # library for plots
library(scales)  # library for transformed axis in plots
library(lemon) # library for facet plots
library(ggh4x)
library(grid) ; library(gridExtra) # libraries for multi-plots
library(ggpubr)   # Library to annotate plots
library(cowplot)  # library to merge multiple plots with labeling (for paper)

library(DescTools) # library for the calculation of AUC
# ------------------ #



# --- FUNCTIONS --- ####
# Functions for plots
exp_bold <- function(lab) {
  new_lab <- log10(lab)
  do.call(
    expression,
    lapply(paste(new_lab), function(x) bquote(bold("10"^.(x))))
  )
}
Plot_colors_Groups <- function(){
  # Colors for the visualization (colors codes can be found in the file 'Recap Info NHPs VAC2105.xlsx')
  colors <-  c(rgb(128,128,128,maxColorValue = 255),rgb(181,219,69,maxColorValue = 255),rgb(69,181,219,maxColorValue = 255),
               rgb(219,69,181,maxColorValue = 255),rgb(219,107,69,maxColorValue = 255),rgb(80,47,168,maxColorValue = 255))
  # names(colors) <- c("Naive","CD40.RBDv Naive","Convalescent","CD40.RBDv Conv","mRNA Conv","CD40.PanCoV Conv")
  names(colors) <- c("Naive","CD40.RBDv - Naive","Convalescent","CD40.RBDv - Conv","mRNA - Conv","CD40.PanCoV - Conv")
  return(colors)
}
Plot_function_Mean_dynamics <- function(data,colors,loq_thresh){
  data$Time <- data$DayPostExpo
  
  # Plot comparing medians of all groups
  Plot_all_group_median <- ggplot(data=data) + 
    geom_rect(data=data.frame(xmin=-1,xmax=30,ymin=0,ymax=loq_thresh),
              aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="gray",alpha=0.3) +
    geom_vline(xintercept = 0,color="red",linetype="longdash",size=0.75) +
    
    geom_line(aes(x=Time,y=Median,color=as.factor(Group)),linewidth=1.5) + 
    scale_color_manual(name="",breaks=names(colors),values=colors) + 
    scale_x_continuous(breaks=seq(0,30,by=5),expand=c(0,0.5)) +
    scale_y_continuous(trans="log10")+
    
    xlab("Time post-exposure (days)") + 
    theme(panel.background = element_rect(fill="white"),
          axis.line = element_line(color="black",size=1.2),
          axis.text =  element_text(color="black",size=11,face="bold"),
          axis.ticks = element_line(color="black",size=1),
          axis.ticks.length = unit(0.1,units = "cm"),
          axis.title = element_text(color="black",size=12,face="bold"),
          legend.key = element_rect(color="white",fill="white",size=0.5),
          strip.background = element_rect(fill="white",color="white"),
          legend.position = "top",
          legend.direction = "horizontal",
          
          legend.text = element_text(face = "bold",color="black",size = 10),
          legend.key.height = unit(0.5,"cm"),
          legend.key.width = unit(1,"cm"),
          strip.text = element_text(size=9,color="black",face="bold")) + 
    guides(color=guide_legend(ncol=2,byrow=TRUE))
  
  return(Plot_all_group_median)
}
Plot_function_Humoral_data <- function(data,colors,pos_thres){
  
  Groups <- unique(data$Group)
  List_plot_groups <- list()
  for(g in 1:length(Groups)){
    # g <- 1
    
    data_group <- data[which(data$Group == Groups[g]),]
    Nb_animal_group <- length(unique(data_group$SubjectID))
    pch_animalID <- rep(c(21,22,23,24,25),ceiling(Nb_animal_group/5))[1:Nb_animal_group]
    
    
    plot_group <- ggplot(data=data_group) + 
      geom_rect(data=data.frame(xmin=-10,xmax=60,ymin=1e-3,ymax=pos_thres),
                aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="gray",alpha=0.3) +
      geom_vline(aes(xintercept = 0,linetype="vaccination"),color="gray40",size=0.5,key_glyph = "path") +
      geom_vline(aes(xintercept = 27,linetype="Infection"),color="red",size=0.5,key_glyph = "path") +
      
      geom_line(aes(x=Time,y=Value,group=as.factor(SubjectID),color=as.factor(Group)),show.legend = FALSE) +
      geom_point(aes(x=Time,y=Value,group=as.factor(SubjectID),fill=as.factor(Group),shape=as.factor(SubjectID)),color="black",cex=2,show.legend = FALSE) +
      facet_grid(.~Group) +
      
      scale_shape_manual(breaks=unique(data_group$SubjectID),values=pch_animalID) +
      scale_color_manual(breaks=unique(data_group$Group),values=colors[Groups[g]]) +
      scale_fill_manual(breaks=unique(data_group$Group),values=colors[Groups[g]]) +
      scale_linetype_manual(name="",breaks = c("vaccination","Infection"),values=c("dotted","longdash")) +
      # scale_linetype_manual(name="",breaks = c("Infection"),values=c("longdash")) + 
      
      
      theme(panel.background = element_rect(fill="white"),
            axis.line = element_line(color="black",size=1.0),
            axis.text =  element_text(color="black",size=9,face="bold"),
            axis.ticks = element_line(color="black",size=0.75),
            axis.ticks.length = unit(0.1,units = "cm"),
            axis.title = element_blank(),
            legend.key = element_rect(color="white",fill="white",size=0.5),
            legend.key.width = unit(units = "cm",x = 2.0),
            legend.key.height = unit(units = "cm",x = 0.25),
            strip.background = element_rect(fill="white",color="white"),
            strip.text = element_text(size=9,color=colors[Groups[g]],face="bold"),
            legend.position = "top") + 
      guides(linetype=guide_legend(override.aes = list(color=c("black","red"),linewidth=1.0)))
    # guides(linetype=guide_legend(override.aes = list(color=c("red"),linewidth=1.0)))
    
    
    List_plot_groups[[Groups[g]]] <- plot_group + theme(axis.title.y = element_blank())
  }
  return(List_plot_groups)
}

# Functions for statistical tests 
TwoGroups_TestComparison <- function(x,y,test_name,...){
  if(test_name == "ttest"){
    # for a t-test we first verify the equality of variance 
    variance_test <- var.test(x = x,y = y)
    if(!is.na(variance_test$p.value)){
      test <- t.test(x = x,y = y ,var.equal = variance_test$p.value > 0.05)
      pvalue <- test$p.value
    }else{
      pvalue <- NA
    }
  }else if(test_name == "wilcoxon"){
    test <- wilcox.test(x=x,y = y)
    pvalue <- test$p.value
  }else{
    pvalue <- NA
  }
  return(pvalue)
}
Global_TestComparison <- function(data,formula,test_name){
  if(test_name == "anova"){
    test <- aov(data=data,formula=formula)
    pvalue <- summary(test)[[1]][["Pr(>F)"]][[1]]
  }else if(test_name == "ks"){
    test <- kruskal.test(data=data,formula)
    pvalue <- test$p.value
  }else{
    pvalue <- NA
  }
  return(pvalue)
}
# ---------------- #








# --- General variables ----
Project_Folder <- "EarlyViralandAntibodyDynPostInfection"
Data_Folder <- "Simulated_data"

Figure_Folder <- paste(Project_Folder,Data_Folder,"Figures",sep="/")
# dir.create(Figure_Folder,recursive = TRUE)   # Uncomment to save figures 

Plot_colors <- Plot_colors_Groups()

Positivity_thresh <- c(BAb = 200, NAb = 0.781250)
LOD <- c(BAb = NA, NAb = 0.1)
# ---------------- #


# --- Download of dataset ----
data_file <- "Simulated_dataset_VL_BAb_NAb_Delta_PostExpo.txt"

Observed_data <- read.table(file = paste(Project_Folder,Data_Folder,data_file,sep="/"),header=TRUE,sep="\t",dec=".")
# Modification of the group name for Convalescent NHPs
Observed_data$Group[which(Observed_data$Group == "Conv")] <- "Convalescent"

# Extraction of Humoral data
Humoral_data <- subset(Observed_data,MeasureType %notin% c("gRNA","sgRNA"))
Humoral_data <- Humoral_data[order(Humoral_data$SubjectID,Humoral_data$MeasureType,Humoral_data$DayPostExpo),]
# Rename of a few columns
colnames(Humoral_data) <- c("SubjectID","Group","StatutImm","DayPostExpo","Value","Censored","ValueCens","IgGType","SampleType",
                            "log10Value","log10ValueCens","obsid","Inoc","Inoc_Voc","Weight","BAb0","NAb0")
# Simplification of the label of IgGType
Humoral_data$IgGType <- with(Humoral_data,ifelse(IgGType == "ACE2 CoV-2 RBD L452R T478K (Delta)","Neutralizing","Binding"))

# Homogenisation of timepoints for Naive group
Humoral_data$DayPostExpo <- with(Humoral_data,case_when(
  Group == "Naive" & DayPostExpo == 23 ~ 22,
  Group == "Naive" & DayPostExpo == 14 ~ 13,
  TRUE ~ DayPostExpo
))

Humoral_data$Group <- factor(Humoral_data$Group,levels = c("Naive","CD40.RBDv - Naive","Convalescent","CD40.RBDv - Conv","mRNA - Conv","CD40.PanCoV - Conv"))
# ---------------- #





# .................................... ----
# --- Plot of individual dynamics ----
# .................................... ----

## > Binding Vs Neutralizing ####
Humoral_data_PostExpo <- Humoral_data[,colnames(Humoral_data) %in% c("SubjectID","Group","DayPostExpo","IgGType","ValueCens")]

Binding_VS_Neutralizing_data <- reshape2::dcast(data=Humoral_data_PostExpo,formula = SubjectID + Group + DayPostExpo ~ IgGType,value.var = "ValueCens")

Colors <- Plot_colors_Groups()

Plot_Binding_Vs_Neutralizing <- ggplot(data=Binding_VS_Neutralizing_data) + 
  geom_point(aes(x=Binding,y=Neutralizing,color=as.factor(Group)),cex=2.0) + 
  
  scale_color_manual(name="Groups",breaks = names(Colors),values=Colors,labels=names(Colors),drop=FALSE) + 
  scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x,n=6),labels = exp_bold) + xlab("anti-RBD binding IgG (AU/mL)") +
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=6),labels = exp_bold) + ylab("ACE2/RBD binding inhibit° (AU/mL)") +
  coord_cartesian(xlim=c(1E1,1.2E6),ylim=c(0.5E-1,1E3),expand = 0) + 
  theme_bw() + 
  guides(color=guide_legend(override.aes = list(linewidth=2))) +
  
  theme(axis.line = element_line(color="black",size=1.0),
        axis.text =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0.1,units = "cm"),
        axis.title = element_text(color="black",size=12,face="bold"),
        legend.background = element_rect(color="black",linewidth=1.0),
        legend.key = element_rect(color="white",fill="white",size=0.5),
        legend.position.inside = c(0.25,0.75),legend.position = "inside",
        legend.title = element_text(face = "bold"),
        legend.text = element_text(face="bold"),
        legend.key.height = unit(units = "cm",x = 0.5),
        panel.grid.major = element_line(colour = "gray90"),
        panel.grid.minor = element_blank())
# ggexport(Plot_Binding_Vs_Neutralizing,filename = paste(Figure_Folder,"Binding_Vs_Neutralizing_PAPER.png",sep="/"),height = 4*300,width=5*300,res = 300)



## > Binding IgG ####
AntiRBD_IgG_data <- subset(Humoral_data,IgGType == "Binding")
AntiRBD_IgG_data$Time <- AntiRBD_IgG_data$DayPostExpo

List_plots_antiRBD_IgG <- Plot_function_Humoral_data(data=AntiRBD_IgG_data,colors=Plot_colors,pos_thres=Positivity_thresh["BAb"])
List_plots_antiRBD_IgG <- lapply(List_plots_antiRBD_IgG,function(plot) plot + 
                                   geom_hline(yintercept = Positivity_thresh["BAb"],color="black",linetype="dashed") +
                                   scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=8),labels = exp_bold) + 
                                   scale_x_continuous(breaks=seq(0,30,by=5)) +
                                   theme(legend.position = "none") +
                                   coord_cartesian(ylim=c(0.085,1e6),xlim=c(-1,23),expand = 0)) 

Plot_Yaxis <- textGrob(expression(bold('Anti-RBD binding IgG (AU/mL)')),rot = 90, gp = gpar(col = "black", fontsize = 10))
Plot_Xaxis <- textGrob(expression(bold('Time post-exposure (days)')),rot = 0,gp = gpar(col = "black", fontsize = 10))

PlotGrid_antiRBD_IgG <- annotate_figure(p=ggarrange(plotlist = List_plots_antiRBD_IgG,nrow=1),left = Plot_Yaxis,bottom = Plot_Xaxis)
# ggexport(PlotGrid_antiRBD_IgG,filename = paste(Figure_Folder,"ACE2RBD_Delta_dynamics_PostExpo.png",sep="/"),height = 3*300,width=15*300,res = 300)


## > Neutralizing IgG ####
ACE2RBD_Inhib_data <- subset(Humoral_data,IgGType == "Neutralizing")
ACE2RBD_Inhib_data$Time <- ACE2RBD_Inhib_data$DayPostExpo

List_plots_ACE2RBD <- Plot_function_Humoral_data(data=ACE2RBD_Inhib_data,colors=Plot_colors,pos_thres=Positivity_thresh["NAb"])
List_plots_ACE2RBD <- lapply(List_plots_ACE2RBD,function(plot) plot + 
                               geom_hline(yintercept = Positivity_thresh["NAb"],color="black",linetype="dashed") +
                               scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=8),labels = exp_bold) + 
                               scale_x_continuous(breaks=seq(0,30,by=5)) +
                               theme(legend.position = "none") +
                               coord_cartesian(ylim=c(0.085,1e6),xlim=c(-1,23),expand = 0))

Plot_Yaxis <- textGrob(expression(bold('ACE2/RBD bind. inh. (AU/mL)')),rot = 90, gp = gpar(col = "black", fontsize = 10))
Plot_Xaxis <- textGrob(expression(bold('Time post-exposure (days)')),rot = 0,gp = gpar(col = "black", fontsize = 10))

PlotGrid_ACE2RBD <- annotate_figure(p=ggarrange(plotlist = List_plots_ACE2RBD,nrow=1),left = Plot_Yaxis,bottom = Plot_Xaxis)
# ggexport(PlotGrid_ACE2RBD,filename = paste(Figure_Folder,"ACE2RBD_Delta_dynamics_PostExpo.png",sep="/"),height = 3*300,width=15*300,res = 300)



## > Merged Plot ####
Merged_IndFigure <- plot_grid(PlotGrid_antiRBD_IgG,PlotGrid_ACE2RBD,
                              labels=c("A","B"),label_fontface = "bold",ncol=1)
# ggsave(Merged_IndFigure,filename = paste(Figure_Folder,"Individual_Ab_dynamics_Merged_Figures_PAPER.png",sep="/"),height = 6,width=15)
# ---------------- #




# .................................... ----
# --- Plot of the mean dynamics ----
# .................................... ----

# Estimation of the median dynamics 
Summary_Antibody_dynamics <- ddply(.data=Humoral_data,.variables = .(IgGType,Group,DayPostExpo),summarise,
                                   Mean=mean(Value,na.rm=TRUE),
                                   Median=median(Value,na.rm=TRUE),
                                   Min = min(Value,na.rm=TRUE), 
                                   Max = max(Value,na.rm=TRUE))



## > Binding IgG ####
Plot_antiRBD_BindingIgG <- Plot_function_Mean_dynamics(data = subset(Summary_Antibody_dynamics,IgGType == "Binding"),
                                                                   colors=Plot_colors,loq_thresh = Positivity_thresh["BAb"]) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n = 8),labels = exp_bold) + 
  geom_hline(yintercept = Positivity_thresh["BAb"],color="black",linetype="dashed") + 
  ylab("Anti-RBD Binding IgG (AU/mL)") +
  scale_x_continuous(breaks=seq(0,30,by=5)) + xlab("Time post-exposure (days)") +
  coord_cartesian(ylim=c(0.085,2e6),xlim=c(-0.25,22),expand = 0) + 
  theme(legend.position.inside = c(0.65,0.20),
        plot.title = element_text(face = "bold",color="black",size=13,hjust = 0)) 
# ggsave(Plot_antiRBD_BindingIgG,filename = paste(Figure_Folder,"Median_Binding_Ab_dynamics_PAPER.png",sep="/"),dpi=300,height = 3.0,width=6)



## > Neutralizing IgG ####
Plot_ACE2RBD_binding_inhibition <- Plot_function_Mean_dynamics(data = subset(Summary_Antibody_dynamics,IgGType == "Neutralizing"),
                                                                           colors=Plot_colors,loq_thresh =  Positivity_thresh["NAb"]) + 
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x,n=8),labels = exp_bold) + 
  geom_hline(yintercept =  Positivity_thresh["NAb"],color="black",linetype="dashed") + 
  ylab("ACE2/RBD bind. inh. (AU/mL)") +
  scale_x_continuous(breaks=seq(0,30,by=5)) + xlab("Time post-exposure (days)") +
  coord_cartesian(ylim=c(0.085,2e6),xlim=c(-0.25,22),expand = 0) + 
  theme(legend.position.inside  = c(0.65,0.80),
        plot.title = element_text(face = "bold",color="black",size=13,hjust = 0)) 
# ggsave(Plot_ACE2RBD_binding_inhibition,filename = paste(Figure_Folder,"Median_Neutralizing_Ab_dynamics_PAPER.png",sep="/"),dpi=300,height = 3,width=6)



## > Merged plot ####
Merged_Figure <- plot_grid(Plot_antiRBD_BindingIgG,Plot_ACE2RBD_binding_inhibition,
                           labels=c("E","F"),label_fontface = "bold",ncol=2)
# ggsave(Merged_Figure,filename = paste(Figure_Folder,"Median_Ab_dynamics_Merged_Figures_PAPER.png",sep="/"),height = 3.5,width=12, dpi=300)

# ---------------- #





# .................................... ----
# --- Calculation of Ab dynamics characteristics (by group) ----
# .................................... ----
# For statistical tests performed on VL descriptors:
#  - For AUC and value at baseline, normal distribution is considered (verified by Shapiro-Wilk normality test)
#        -> perform T-test for 2-2 mean comparisons (welch's or student t-test according to variance equality)
#        -> perform ANOVA for global mean comparisons


Humoral_data_PostExpo <- subset(Humoral_data,DayPostExpo >=0)
# Re-order dataset by time points to rightly calculate AUC
Humoral_data_PostExpo <- Humoral_data_PostExpo[order(Humoral_data_PostExpo$DayPostExpo),]

Overall_Ab_descriptors <- data.frame(IgGType=character(),
                                     Group=character(),SubjectID=character(),
                                     Descriptor=character(),Value=numeric())


## 1. Positivity at Exposure ----
# Estimation of the number of NHP with values above positivity threshold at exposure
Humoral_data_PostExpo$Positivity <- with(Humoral_data_PostExpo,
                                         case_when(IgGType == "Binding" & Value >= Positivity_thresh["BAb"] ~ 1,
                                                   IgGType == "Neutralizing" & Value >= Positivity_thresh["NAb"] ~ 1,
                                                   TRUE ~ 0))
Distribution_Postivity_Exposure <- ddply(.data=subset(Humoral_data_PostExpo,DayPostExpo == 0),.variables = .(IgGType,Group),summarise,
                                         N = length(Positivity),
                                         Npos = sum(Positivity,na.rm=TRUE),
                                         PctPos = sum(Positivity,na.rm=TRUE)*100/length(Positivity))




## 2. AUC of Ab dynamics ----

### > Calculation of descriptors ----
# For Binding antibodies, we calculated AUC on log10 value (all values are higher than 1 AU/mL)
# For Neutralizing antibodies, we calculated AUC above -1 log10 AU/mL  (i.e. calculated AUC using yij + 1) 

AUC_Ab_Binding <- ddply(.data=subset(Humoral_data_PostExpo,DayPostExpo <= 23 & IgGType == "Binding"),
                        .variables = .(IgGType,SubjectID,Group),summarize,
                        AUClog10 = AUC(x = DayPostExpo,y=log10(Value),na.rm = TRUE),
                        AUC = AUC(x = DayPostExpo,y=Value,na.rm = TRUE))

AUC_Ab_Neutralizing <- ddply(.data=subset(Humoral_data_PostExpo,DayPostExpo <= 23 & IgGType == "Neutralizing"),
                             .variables = .(IgGType,SubjectID,Group),summarize,
                             AUClog10 = AUC(x = DayPostExpo,y=log10(Value)+1,na.rm = TRUE),
                             AUC = AUC(x = DayPostExpo,y=Value,na.rm = TRUE))

AUC_Ab <- rbind(AUC_Ab_Binding,AUC_Ab_Neutralizing)

Overall_Ab_descriptors <- rbind(Overall_Ab_descriptors,melt(data=AUC_Ab,id.vars = c("IgGType","Group","SubjectID"),variable.name = "Descriptor",value.name = "Value"))

Distribution_AUC <- ddply(.data = AUC_Ab,.variables = .(IgGType,Group),summarize,
                          N = length(SubjectID),
                          Mean_AUC = mean(AUC,na.rm=TRUE),
                          ICMIN_AUC = quantile(AUC,na.rm=TRUE,probs=c(0.025)),
                          ICMAX_AUC = quantile(AUC,na.rm=TRUE,probs=c(0.975)))


### > Global Statistical test Comparison ----
Global_TestComparison_AUC <- expand.grid(IgGType=c("Binding","Neutralizing"),Descriptor=c("AUC","AUClog10"))
Global_TestComparison_AUC$Test <- "anova"
Global_TestComparison_AUC$Pvalue <- sapply(1:nrow(Global_TestComparison_AUC),function(i) Global_TestComparison(data=subset(AUC_Ab,IgGType == Global_TestComparison_AUC$IgGType[i]),
                                                                                                               formula = formula(ifelse(Global_TestComparison_AUC$Descriptor[i] == "AUC",'AUC ~ Group','AUClog10~Group')),
                                                                                                               test_name = Global_TestComparison_AUC$Test[i]))
Global_TestComparison_AUC$Adj.Pvalue <- p.adjust(Global_TestComparison_AUC$Pvalue,method="BH")


### > Two-sample comparison ----
TwoSample_test_AUC <- data.frame(gtools::combinations(n = 6, r = 2, repeats.allowed = F, v = as.character(unique(AUC_Ab$Group))))
TwoSample_test_AUC <- rbind(data.frame(IgGType = "Binding",Group1 = TwoSample_test_AUC$X1, Group2 = TwoSample_test_AUC$X2),
                            data.frame(IgGType = "Neutralizing",Group1 = TwoSample_test_AUC$X1, Group2 = TwoSample_test_AUC$X2))
TwoSample_test_AUC$Pvalue.AUC <- NA
TwoSample_test_AUC$Pvalue.AUClog10 <- NA
for(i in 1:nrow(TwoSample_test_AUC)){
  # i <- 1
  
  tmp_data_G1 <- subset(AUC_Ab,IgGType == TwoSample_test_AUC$IgGType[i] & Group == TwoSample_test_AUC$Group1[i] & !is.na(AUC))
  tmp_data_G2 <- subset(AUC_Ab,IgGType == TwoSample_test_AUC$IgGType[i] & Group == TwoSample_test_AUC$Group2[i] & !is.na(AUC))
  
  if(nrow(tmp_data_G1) > 1 & nrow(tmp_data_G2) > 1){
    TwoSample_test_AUC$Pvalue.AUC[i] <- TwoGroups_TestComparison(x=tmp_data_G1$AUC,y=tmp_data_G2$AUC,test_name = "ttest")
    TwoSample_test_AUC$Pvalue.AUClog10[i] <- TwoGroups_TestComparison(x=tmp_data_G1$AUClog10,y=tmp_data_G2$AUClog10,test_name = "ttest")
  }
}
TwoSample_test_AUC$Adj.Pvalue.AUC <- NA ; TwoSample_test_AUC$Adj.Pvalue.AUClog10 <- NA
TwoSample_test_AUC$Adj.Pvalue.AUC[which(TwoSample_test_AUC$IgGType == "Binding")] <- p.adjust(TwoSample_test_AUC$Pvalue.AUC[which(TwoSample_test_AUC$IgGType == "Binding")],method="BH")
TwoSample_test_AUC$Adj.Pvalue.AUC[which(TwoSample_test_AUC$IgGType == "Neutralizing")] <- p.adjust(TwoSample_test_AUC$Pvalue.AUC[which(TwoSample_test_AUC$IgGType == "Neutralizing")],method="BH")
TwoSample_test_AUC$Significant.AUC <- 1*(TwoSample_test_AUC$Adj.Pvalue.AUC <= 0.05)

TwoSample_test_AUC$Adj.Pvalue.AUClog10[which(TwoSample_test_AUC$IgGType == "Binding")] <- p.adjust(TwoSample_test_AUC$Pvalue.AUClog10[which(TwoSample_test_AUC$IgGType == "Binding")],method="BH")
TwoSample_test_AUC$Adj.Pvalue.AUClog10[which(TwoSample_test_AUC$IgGType == "Neutralizing")] <- p.adjust(TwoSample_test_AUC$Pvalue.AUClog10[which(TwoSample_test_AUC$IgGType == "Neutralizing")],method="BH")
TwoSample_test_AUC$Significant.AUClog10 <- 1*(TwoSample_test_AUC$Adj.Pvalue.AUClog10 <= 0.05)




## 3. Value at Exposure ----
### > Calculation of descriptors ----
Baseline_Ab <- ddply(.data=subset(Humoral_data_PostExpo,DayPostExpo == 0),
                     .variables = .(IgGType,SubjectID,Group),summarize,
                     Ab0 = Value,Ab0log10 = log10(Value))

Overall_Ab_descriptors <- rbind(Overall_Ab_descriptors,melt(data=Baseline_Ab,id.vars = c("IgGType","Group","SubjectID"),variable.name = "Descriptor",value.name = "Value"))

Distribution_Baseline_Ab <- ddply(.data = Baseline_Ab,.variables = .(IgGType,Group),summarize,
                                  N = length(SubjectID),
                                  Mean_Ab0 = mean(Ab0,na.rm=TRUE),
                                  ICMIN_Ab0 = quantile(Ab0,na.rm=TRUE,probs=c(0.025)),
                                  ICMAX_Ab0 = quantile(Ab0,na.rm=TRUE,probs=c(0.975)))

Distribution_Baseline_Ab$Label <- with(Distribution_Baseline_Ab,
                                       paste(round(Mean_Ab0,digits=2)," [",round(ICMIN_Ab0,digits=2),
                                             " ; ",round(ICMAX_Ab0,digits=2),"]",sep=""))

### > Global Statistical test Comparison ----
Global_TestComparison_Ab0 <- expand.grid(IgGType=c("Binding","Neutralizing"),Descriptor=c("Ab0","Ab0log10"))
Global_TestComparison_Ab0$Test <- "anova"
Global_TestComparison_Ab0$Pvalue <- sapply(1:nrow(Global_TestComparison_Ab0),function(i) Global_TestComparison(data=subset(Baseline_Ab,IgGType == Global_TestComparison_AUC$IgGType[i]),
                                                                                                               formula = formula(ifelse(Global_TestComparison_Ab0$Descriptor[i] == "Ab0",'Ab0 ~ Group','Ab0log10~Group')),
                                                                                                               test_name = Global_TestComparison_AUC$Test[i]))
Global_TestComparison_Ab0$Adj.Pvalue <- p.adjust(Global_TestComparison_Ab0$Pvalue,method="BH")


### > Two-sample comparison ----
TwoSample_test_Ab0 <- data.frame(gtools::combinations(n = 6, r = 2, repeats.allowed = F, v = as.character(unique(Baseline_Ab$Group))))
TwoSample_test_Ab0 <- rbind(data.frame(IgGType = "Binding",Group1 = TwoSample_test_Ab0$X1, Group2 = TwoSample_test_Ab0$X2),
                            data.frame(IgGType = "Neutralizing",Group1 = TwoSample_test_Ab0$X1, Group2 = TwoSample_test_Ab0$X2))
TwoSample_test_Ab0$Pvalue.Ab0 <- NA
TwoSample_test_Ab0$Pvalue.Ab0log10 <- NA
for(i in 1:nrow(TwoSample_test_Ab0)){
  # i <- 1
  
  tmp_data_G1 <- subset(Baseline_Ab,IgGType == TwoSample_test_Ab0$IgGType[i] & Group == TwoSample_test_Ab0$Group1[i] & !is.na(Ab0))
  tmp_data_G2 <- subset(Baseline_Ab,IgGType == TwoSample_test_Ab0$IgGType[i] & Group == TwoSample_test_Ab0$Group2[i] & !is.na(Ab0))
  
  if(nrow(tmp_data_G1) > 1 & nrow(tmp_data_G2) > 1){
    TwoSample_test_Ab0$Pvalue.Ab0[i] <- TwoGroups_TestComparison(x=tmp_data_G1$Ab0,y=tmp_data_G2$Ab0,test_name = "ttest")
    TwoSample_test_Ab0$Pvalue.Ab0log10[i] <- TwoGroups_TestComparison(x=tmp_data_G1$Ab0log10,y=tmp_data_G2$Ab0log10,test_name = "ttest")
  }
}
TwoSample_test_Ab0$Adj.Pvalue.Ab0 <- NA ; TwoSample_test_Ab0$Adj.Pvalue.Ab0log10 <- NA
TwoSample_test_Ab0$Adj.Pvalue.Ab0[which(TwoSample_test_Ab0$IgGType == "Binding")] <- p.adjust(TwoSample_test_Ab0$Pvalue.Ab0[which(TwoSample_test_Ab0$IgGType == "Binding")],method="BH")
TwoSample_test_Ab0$Adj.Pvalue.Ab0[which(TwoSample_test_Ab0$IgGType == "Neutralizing")] <- p.adjust(TwoSample_test_Ab0$Pvalue.Ab0[which(TwoSample_test_Ab0$IgGType == "Neutralizing")],method="BH")
TwoSample_test_Ab0$Significant.Ab0 <- 1*(TwoSample_test_Ab0$Adj.Pvalue.Ab0 <= 0.05)

TwoSample_test_Ab0$Adj.Pvalue.Ab0log10[which(TwoSample_test_Ab0$IgGType == "Binding")] <- p.adjust(TwoSample_test_Ab0$Pvalue.Ab0log10[which(TwoSample_test_Ab0$IgGType == "Binding")],method="BH")
TwoSample_test_Ab0$Adj.Pvalue.Ab0log10[which(TwoSample_test_Ab0$IgGType == "Neutralizing")] <- p.adjust(TwoSample_test_Ab0$Pvalue.Ab0log10[which(TwoSample_test_Ab0$IgGType == "Neutralizing")],method="BH")
TwoSample_test_Ab0$Significant.Ab0log10 <- 1*(TwoSample_test_Ab0$Adj.Pvalue.Ab0log10 <= 0.05)






## 4. Plot of results ----
Plot_colors <- Plot_colors_Groups()
Overall_Ab_descriptors$Group <- factor(Overall_Ab_descriptors$Group,levels = names(Plot_colors))
Overall_Ab_descriptors$Descriptor <- as.character(Overall_Ab_descriptors$Descriptor)
Overall_Ab_descriptors$Descriptor <- factor(Overall_Ab_descriptors$Descriptor,levels=c("Ab0","Ab0log10","AUC","AUClog10"))

# Calculation of the distribution of each descriptor
Distribution_Descriptors <- ddply(.data = Overall_Ab_descriptors,.variables = .(IgGType,Group,Descriptor),summarize,
                                  Mean = mean(Value,na.rm=TRUE),
                                  Median = median(Value,na.rm=TRUE),
                                  ICMIN=quantile(Value,0.025,na.rm=TRUE),
                                  IC25=quantile(Value,0.25,na.rm=TRUE),
                                  IC75=quantile(Value,0.75,na.rm=TRUE),
                                  ICMAX=quantile(Value,0.975,na.rm=TRUE))


Global_TestComparison <- rbind(Global_TestComparison_AUC,Global_TestComparison_Ab0)
Global_TestComparison$Significant <- with(Global_TestComparison,1*(Adj.Pvalue<=0.05))
Global_TestComparison$Label <- with(Global_TestComparison,
                                    case_when(Test == "anova" & Adj.Pvalue < 1e-3 ~ "ANOVA pvalue < 0.001",
                                              Test == "ks" & Adj.Pvalue < 1e-3 ~ "KS pvalue < 0.001",
                                              Test == "anova" & Adj.Pvalue >= 1e-3 ~ paste("ANOVA pvalue = ",round(Adj.Pvalue,digits=3),sep=""),
                                              Test == "ks" & Adj.Pvalue >= 1e-3 ~ paste("KS pvalue = ",round(Adj.Pvalue,digits=3),sep="")))
Global_TestComparison$x <- with(Global_TestComparison,
                                case_when(Descriptor == "Ab0log10" & IgGType == "Binding" ~ 1,
                                          Descriptor == "Ab0log10" & IgGType == "Neutralizing" ~ -1,
                                          Descriptor == "AUClog10" & IgGType == "Binding" ~ 60,
                                          Descriptor == "AUClog10" & IgGType == "Neutralizing" ~ 20))
Global_TestComparison$y <- 7.2

Labels_Descriptors <- c("Value at exposure \n [in log10 AU/mL]",
                        "Area under immune curve \n [in log10 AU.day/mL]")
names(Labels_Descriptors) <- c("Ab0log10","AUClog10")


### > Binding IgG  ----
GeomPointRange_Binding <- ggplot(data=subset(Distribution_Descriptors,IgGType == "Binding" & Descriptor %in% c("Ab0log10","AUClog10"))) + 
  geom_pointrange(aes(x=Mean,xmin=ICMIN,xmax=ICMAX,y=Group,color=Group),shape=21,stroke=1.5,cex=0.60,fill="white",linewidth=1.0,show.legend = FALSE) + 
  facet_rep_grid(.~Descriptor,scales="free",repeat.tick.labels=FALSE,labeller=labeller(Descriptor = Labels_Descriptors)) + 
  geom_text(data=subset(Global_TestComparison,IgGType == "Binding" & Descriptor %in% c("Ab0log10","AUClog10")),
            aes(x=x,y=y,label=Label),hjust = 0,size = 3.5,fontface = "plain") +
  scale_color_manual(values = Plot_colors) + 
  
  facetted_pos_scales(x= list(
    Descriptor == "Ab0log10" ~ scale_x_continuous(breaks =seq(0,6,by=1),limits=c(1,6)),
    Descriptor == "AUClog10" ~ scale_x_continuous(breaks = seq(0,120,by=20), limits=c(60,120)))) + 
  
  coord_cartesian(ylim=c(1,7)) +
  ggtitle("anti-RBD binding IgG") + 
  theme_bw() + 
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(color="black",size=0.75),
        axis.text.x =  element_text(color="black",size=10),
        axis.text.y =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0,units = "cm"),
        axis.title = element_blank(),
        plot.title = element_text(color="black",face = "bold.italic",size = 12),
        strip.background = element_rect(fill=alpha("darkblue",alpha = 0.2),color="darkblue",linewidth = 0.75),
        strip.text = element_text(color="black",face="bold",size=11,margin = margin(b=1,t=1)),
        legend.position = "none") 


### > Neutralizing IgG  ---- 
GeomPointRange_Neutralizing <- ggplot(data=subset(Distribution_Descriptors,IgGType == "Neutralizing" & Descriptor %in% c("Ab0log10","AUClog10"))) + 
  geom_pointrange(aes(x=Mean,xmin=ICMIN,xmax=ICMAX,y=Group,color=Group),shape=21,stroke=1.5,cex=0.60,fill="white",linewidth=1.0,show.legend = FALSE) + 
  facet_rep_grid(.~Descriptor,scales="free",repeat.tick.labels=FALSE,labeller=labeller(Descriptor = Labels_Descriptors)) + 
  geom_text(data=subset(Global_TestComparison,IgGType == "Neutralizing" & Descriptor %in% c("Ab0log10","AUClog10")),
            aes(x=x,y=y,label=Label),hjust = 0,size = 3.5,fontface = "plain") +
  scale_color_manual(values = Plot_colors) + 
  
  facetted_pos_scales(x= list(
    Descriptor == "Ab0log10" ~ scale_x_continuous(breaks =seq(-1,6,by=1),limits=c(-1,3)),
    Descriptor == "AUClog10" ~ scale_x_continuous(breaks = seq(0,120,by=20), limits=c(20,80)))) + 
  
  coord_cartesian(ylim=c(1,7)) +
  ggtitle("ACE2/RBD binding inhibition") + 
  theme_bw() + 
  theme(panel.background = element_rect(fill="white"),
        axis.line = element_line(color="black",size=0.75),
        axis.text.x =  element_text(color="black",size=10),
        axis.text.y =  element_text(color="black",size=10,face="bold"),
        axis.ticks = element_line(color="black",size=1),
        axis.ticks.length = unit(0,units = "cm"),
        axis.title = element_blank(),
        plot.title = element_text(color="black",face = "bold.italic",size = 12),
        strip.background = element_rect(fill=alpha("darkblue",alpha = 0.2),color="darkblue",linewidth = 0.75),
        strip.text = element_text(color="black",face="bold",size=11,margin = margin(b=1,t=1)),
        legend.position = "none") 


### > Merged Plot ---- 
GeomPointRange_Antibodies <- plot_grid(GeomPointRange_Binding,NULL,GeomPointRange_Neutralizing,nrow=1,align="v",labels=c("C","","D"),rel_widths = c(1,0.1,1))
# ggsave(plot=GeomPointRange_Antibodies,filename = paste(Figure_Folder,"Antibody_Descriptors_PAPER.png",sep ="/"),height=3,width=12,dpi=300)
