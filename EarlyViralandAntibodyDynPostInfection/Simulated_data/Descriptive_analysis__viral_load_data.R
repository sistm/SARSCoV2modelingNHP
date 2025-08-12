# --------------- # 
# DESCRIPTION: Descriptive analysis of viral load dynamics used in the modeling work
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
library(grid) ; library(gridExtra) # libraries for multi-plots
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
  names(colors) <- c("Naive","CD40.RBDv - Naive","Convalescent","CD40.RBDv - Conv","mRNA - Conv","CD40.PanCoV - Conv")
  return(colors)
}
Plot_function_Viral_Load <- function(data,colors,loq_thresh){
  
  Groups <- unique(data$Group)
  List_plot_groups <- list()
  for(g in 1:length(Groups)){
    
    # g <- 1
    data_group <- subset(data,Group == Groups[g])
    Nb_animal_group <- length(unique(data_group$SubjectID))
    pch_SubjectID <- rep(c(21,22,23,24,25),ceiling(Nb_animal_group/5))[1:Nb_animal_group]
    
    plot_group <- ggplot(data=data_group) + 
      geom_rect(data=data.frame(xmin=-0.5,xmax=max(data$Time)+0.5,ymin=0.8*min(data$log10VLcens),ymax=log10(loq_thresh)),
                aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="gray",alpha=0.3) +
      geom_vline(xintercept = 0,color="red",linetype="longdash",linewidth=0.5) +
      
      geom_line(aes(x=Time,y=log10VLcens,group=as.factor(SubjectID),color=as.factor(Group))) + 
      geom_point(aes(x=Time,y=log10VLcens,group=as.factor(SubjectID),fill=as.factor(Group),shape=as.factor(SubjectID)),color="black",cex=2) + 
      facet_grid(.~Group) + 
      
      scale_shape_manual(breaks=unique(data_group$SubjectID),values=pch_SubjectID) + 
      scale_color_manual(breaks=unique(data_group$Group),values=colors[Groups[g]]) + 
      scale_fill_manual(breaks=unique(data_group$Group),values=colors[Groups[g]]) + 
      scale_x_continuous(breaks=seq(0,25,by=5),expand=c(0,0.5)) +
      scale_y_continuous(breaks=seq(0,10,by=1)) +
      coord_cartesian(xlim=c(0,NA)) +
      xlab("Time post-exposure (days)") + 
      theme(panel.background = element_rect(fill="white"),
            axis.line = element_line(color="black",linewidth=1.0),
            axis.text =  element_text(color="black",size=9,face="bold"),
            axis.ticks = element_line(color="black",linewidth=0.75),
            axis.ticks.length = unit(0.1,units = "cm"),
            axis.title = element_text(color="black",size=10,face="bold"),
            legend.key = element_rect(color="white",fill="white",size=0.5),
            strip.background = element_rect(fill="white",color="white"),
            strip.text = element_text(size=9,color=colors[Groups[g]],face="bold"),
            legend.position = "none")
    
    List_plot_groups[[Groups[g]]] <- plot_group
  }
  return(List_plot_groups)
}
Plot_function_Mean_Viral_Load <- function(data,colors,loq_thresh){
  
  # Plot comparing medians of all groups
  Plot_all_group_median <- ggplot(data=data) + 
    geom_rect(data=data.frame(xmin=-0.5,xmax=max(data$Time)+0.5,ymin=0.8*min(data$Median),ymax=loq_thresh),
              aes(xmin=xmin,xmax=xmax,ymin=ymin,ymax=ymax),fill="gray",alpha=0.3) +
    geom_vline(xintercept = 0,color="red",linetype="longdash",linewidth=0.75) +
    
    geom_line(aes(x=Time,y=Median,color=as.factor(Group)),linewidth=1.5) + 
    scale_color_manual(name="",breaks=names(colors),values=colors) + 
    scale_x_continuous(breaks=seq(0,25,by=5),expand=c(0,0.5)) +
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),labels = exp_bold) +
    
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
          legend.key.width = unit(1,"cm")) + 
    guides(color=guide_legend(ncol=2,byrow=TRUE))
  
  return(Plot_all_group_median)
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
# ---------------- #





# --- Download of dataset ----
data_file <- "Simulated_dataset_VL_BAb_NAb_Delta_PostExpo.txt"

LOD_gRNA <- 476 # in copies/mL
LOD_sgRNA <- 749 # in copies/mL
LOQ_gRNA <- LOD_gRNA*10
LOQ_sgRNA <- LOD_sgRNA*10

Observed_data <- read.table(file = paste(Project_Folder,Data_Folder,data_file,sep="/"),header=TRUE,sep="\t",dec=".")
# Modification of the group name for Convalescent NHPs
Observed_data$Group[which(Observed_data$Group == "Conv")] <- "Convalescent"

# Extraction of viral load data
Viral_Load_data <- subset(Observed_data,MeasureType %in% c("gRNA","sgRNA"))
# Rename of a few columns
colnames(Viral_Load_data) <- c("SubjectID","Group","StatutImm","Time","VL","Censored","VLcens","VLType","SampleType",
                               "log10VL","log10VLcens","obsid","Inoc","Inoc_Voc","Weight","BAb0","NAb0")
Viral_Load_data$Group <- factor(Viral_Load_data$Group,levels = c("Naive","CD40.RBDv - Naive","Convalescent","CD40.RBDv - Conv","mRNA - Conv","CD40.PanCoV - Conv"))
# ---------------- #




# .................................... ----
# --- Plot of individual dynamics ----
# .................................... ----

## > gRNA VL in Trachea ####
Viral_Load_gRNA_Trachea_data <- subset(Viral_Load_data,VLType == "gRNA" & SampleType == "Tracheal fluid" & Time >=0)
List_plots_gRNA_Trachea <- Plot_function_Viral_Load(data=Viral_Load_gRNA_Trachea_data,colors = Plot_colors,loq_thresh=LOQ_gRNA)

List_plots_gRNA_Trachea <- lapply(List_plots_gRNA_Trachea,function(plot) plot + scale_y_continuous(breaks = seq(1,10,by=1)) + 
                                    coord_cartesian(ylim=c(2.5,9.1),xlim=c(-0.5,30),expand = 0) + 
                                    theme(axis.title = element_blank()))
Plot_Yaxis <- textGrob(expression(bold('gRNA (log' [10]* ' cp/mL)')),rot = 90,gp = gpar(col = "black", fontsize = 10))
Plot_Xaxis <- textGrob(expression(bold('Time post-exposure (days)')),rot = 0,gp = gpar(col = "black", fontsize = 10))
Plot_Title <- textGrob(expression(bold('Trachea')),rot = 0,x=0,hjust = 0,gp = gpar(col = "black", fontsize = 12))

PlotGrid_gRNA_Trachea <- grid.arrange(grobs=List_plots_gRNA_Trachea,left=Plot_Yaxis,bottom=Plot_Xaxis,top=Plot_Title,nrow=1)
# ggsave(PlotGrid_gRNA_Trachea,filename = paste(Figure_Folder,"VL_gRNA_Trachea_dynamics.png",sep="/"),height = 3,width=15)


## > gRNA VL in Nasopharynx ####
Viral_Load_gRNA_Naso_data <- subset(Viral_Load_data,VLType == "gRNA" & SampleType == "Nasal fluid" & Time >=0)
List_plots_gRNA_Naso <- Plot_function_Viral_Load(data=Viral_Load_gRNA_Naso_data,colors = Plot_colors,loq_thresh=LOQ_gRNA)

List_plots_gRNA_Naso <- lapply(List_plots_gRNA_Naso,function(plot) plot + scale_y_continuous(breaks = seq(1,10,by=1)) +
                                 coord_cartesian(ylim=c(2.5,9.1),xlim=c(-0.5,30),expand = 0) + 
                                 theme(axis.title = element_blank()))
Plot_Yaxis <- textGrob(expression(bold('gRNA (log' [10]* ' cp/mL)')),rot = 90,gp = gpar(col = "black", fontsize = 10))
Plot_Xaxis <- textGrob(expression(bold('Time post-exposure (days)')),rot = 0,gp = gpar(col = "black", fontsize = 10))
Plot_Title <- textGrob(expression(bold('Nasopharynx')),rot = 0,x=0,hjust = 0,gp = gpar(col = "black", fontsize = 12))

PlotGrid_gRNA_Naso <- grid.arrange(grobs=List_plots_gRNA_Naso,left=Plot_Yaxis,bottom=Plot_Xaxis,top=Plot_Title,nrow=1)
# ggsave(PlotGrid_gRNA_Naso,filename = paste(Figure_Folder,"VL_gRNA_Naso_dynamics.png",sep="/"),height = 3,width=15)


## > sgRNA VL in Trachea #### 
Viral_Load_sgRNA_Trachea_data <-  subset(Viral_Load_data,VLType == "sgRNA" & SampleType == "Tracheal fluid" & Time >=0)
List_plots_sgRNA_Trachea <- Plot_function_Viral_Load(data=Viral_Load_sgRNA_Trachea_data,colors = Plot_colors,loq_thresh=LOQ_sgRNA)

List_plots_sgRNA_Trachea <- lapply(List_plots_sgRNA_Trachea,function(plot) plot + scale_y_continuous(breaks = seq(1,10,by=1)) +
                                     coord_cartesian(ylim=c(2.5,9.1),xlim=c(-0.5,30),expand = 0) + 
                                     theme(axis.title = element_blank()))
Plot_Yaxis <- textGrob(expression(bold('sgRNA (log' [10]* ' cp/mL)')),rot = 90,gp = gpar(col = "black", fontsize = 10))
Plot_Xaxis <- textGrob(expression(bold('Time post-exposure (days)')),rot = 0,gp = gpar(col = "black", fontsize = 10))
Plot_Title <- textGrob(expression(bold('Trachea')),rot = 0,x=0,hjust = 0,gp = gpar(col = "black", fontsize = 12))

PlotGrid_sgRNA_Trachea <- grid.arrange(grobs=List_plots_sgRNA_Trachea,left=Plot_Yaxis,bottom=Plot_Xaxis,top=Plot_Title,nrow=1)
# ggsave(PlotGrid_sgRNA_Trachea,filename = paste(Figure_Folder,"VL_sgRNA_Trachea_dynamics.png",sep="/"),height = 3,width=15)


## > sgRNA VL in Nasopharynx ####
Viral_Load_sgRNA_Naso_data <- subset(Viral_Load_data,VLType == "sgRNA" & SampleType == "Nasal fluid" & Time >=0)
List_plots_sgRNA_Naso <- Plot_function_Viral_Load(data=Viral_Load_sgRNA_Naso_data,colors = Plot_colors,loq_thresh=LOQ_sgRNA)

List_plots_sgRNA_Naso <- lapply(List_plots_sgRNA_Naso,function(plot) plot + scale_y_continuous(breaks = seq(1,10,by=1)) + 
                                  coord_cartesian(ylim=c(2.5,9.1),xlim=c(-0.5,30),expand = 0) + 
                                  theme(axis.title = element_blank()))
Plot_Yaxis <- textGrob(expression(bold('sgRNA (log' [10]* ' cp/mL)')),rot = 90,gp = gpar(col = "black", fontsize = 10))
Plot_Xaxis <- textGrob(expression(bold('Time post-exposure (days)')),rot = 0,gp = gpar(col = "black", fontsize = 10))
Plot_Title <- textGrob(expression(bold('Nasopharynx')),rot = 0,x=0,hjust = 0,gp = gpar(col = "black", fontsize = 12))

PlotGrid_sgRNA_Naso <- grid.arrange(grobs=List_plots_sgRNA_Naso,left=Plot_Yaxis,bottom=Plot_Xaxis,top=Plot_Title,nrow=1)
# ggsave(PlotGrid_sgRNA_Naso,filename = paste(Figure_Folder,"VL_sgRNA_Naso_dynamics.png",sep="/"),height = 3,width=15)



## > Merged plots (for PAPER) ####
Merged_IndFigure <- plot_grid(PlotGrid_gRNA_Trachea,PlotGrid_gRNA_Naso,PlotGrid_sgRNA_Trachea,PlotGrid_sgRNA_Naso,
                              labels=c("A","B","C","D"),
                              label_fontface = "bold",
                              ncol=1)

# ggsave(Merged_IndFigure,filename = paste(Figure_Folder,"Individual_VL_dynamics_Merged_Figures_PAPER.png",sep="/"),height = 12,width=15)
# ---------------- #



# .................................... ----
# --- Plot of the mean dynamics ----
# .................................... ----

# Estimation of the mean, median, 95th CI, sd and number of data at each time point within each group of treatment 
Viral_Load_summary_dynamics <- ddply(.data=Viral_Load_data,.variables = .(VLType,SampleType,Time,Group),summarise,
                                     N=length(unique(SubjectID)),
                                     Mean=mean(VLcens,na.rm=TRUE),
                                     Median=median(VLcens,na.rm=TRUE),
                                     Q1=quantile(VLcens,probs = c(0.25),na.rm=TRUE),      # 25th percentile
                                     Q3=quantile(VLcens,probs = c(0.75),na.rm=TRUE),      # 75th percentile
                                     ICMIN=quantile(VLcens,probs = c(0.025),na.rm=TRUE),  # 2.5th percentile
                                     ICMAX=quantile(VLcens,probs = c(0.975),na.rm=TRUE))  # 97.5th percentile


## > gRNA VL in Trachea ####
Viral_Load_sumDyn_gRNA_Trachea <- subset(Viral_Load_summary_dynamics,VLType == "gRNA" & SampleType == "Tracheal fluid" & Time >=0)
Plot_sumDyn_gRNA_Trachea <- Plot_function_Mean_Viral_Load(data=Viral_Load_sumDyn_gRNA_Trachea,colors=Plot_colors,loq_thresh=LOQ_gRNA)

Plot_sumDyn_gRNA_Trachea <- Plot_sumDyn_gRNA_Trachea + 
  ylab("Median gRNA (copies/mL)") +
  ggtitle("Trachea") + 
  theme(legend.position.inside = c(0.65,0.80),
        plot.title = element_text(face = "bold",color="black",size=13,hjust = 0)) +
  coord_cartesian(ylim=c(0.8*min(Viral_Load_sumDyn_gRNA_Trachea$Median),1E9),xlim=c(-0.25,30),expand = 0)
# ggsave(Plot_sumDyn_gRNA_Trachea,filename = paste(Figure_Folder,"Median_VL_gRNA_Trachea_dynamics_PAPER.png",sep="/"),dpi=300,height = 3,width=6)


## > gRNA VL in Naso ####
Viral_Load_sumDyn_gRNA_Naso <- subset(Viral_Load_summary_dynamics,VLType == "gRNA" & SampleType == "Nasal fluid" & Time >=0)
Plot_sumDyn_gRNA_Naso <- Plot_function_Mean_Viral_Load(data=Viral_Load_sumDyn_gRNA_Naso,colors=Plot_colors,loq_thresh=LOQ_gRNA)

Plot_sumDyn_gRNA_Naso <- Plot_sumDyn_gRNA_Naso + 
  ylab("Median gRNA (copies/mL)") +
  ggtitle("Nasopharynx") + 
  theme(legend.position.inside = c(0.65,0.80),
        plot.title = element_text(face = "bold",color="black",size=13,hjust = 0)) +
  coord_cartesian(ylim=c(0.8*min(Viral_Load_sumDyn_gRNA_Naso$Median),1E9),xlim=c(-0.25,30),expand = 0)
# ggsave(Plot_sumDyn_gRNA_Naso,filename = paste(Figure_Folder,"Median_VL_gRNA_Naso_dynamics_PAPER.png",sep="/"),dpi=300,height = 3,width=6)


## > sgRNA VL in Trachea ####
Viral_Load_sumDyn_sgRNA_Trachea <- subset(Viral_Load_summary_dynamics,VLType == "sgRNA" & SampleType == "Tracheal fluid" & Time >=0)
Plot_sumDyn_sgRNA_Trachea <- Plot_function_Mean_Viral_Load(data=Viral_Load_sumDyn_sgRNA_Trachea,colors=Plot_colors,loq_thresh=LOQ_sgRNA) 

Plot_sumDyn_sgRNA_Trachea <- Plot_sumDyn_sgRNA_Trachea + 
  ylab("Median sgRNA (copies/mL)") +
  ggtitle("Trachea") + 
  theme(legend.position.inside = c(0.65,0.80),
        plot.title = element_text(face = "bold",color="black",size=13,hjust = 0)) +
  coord_cartesian(ylim=c(0.8*min(Viral_Load_sumDyn_sgRNA_Trachea$Median),1E9),xlim=c(-0.25,30),expand = 0)
# ggsave(Plot_sumDyn_sgRNA_Trachea,filename = paste(Figure_Folder,"Median_VL_sgRNA_Trachea_dynamics_PAPER.png",sep="/"),dpi=300, height = 3,width=6)


## > gRNA VL in Naso ####
Viral_Load_sumDyn_sgRNA_Naso <- subset(Viral_Load_summary_dynamics,VLType == "sgRNA" & SampleType == "Nasal fluid" & Time >=0)
Plot_sumDyn_sgRNA_Naso <- Plot_function_Mean_Viral_Load(data=Viral_Load_sumDyn_sgRNA_Naso,colors=Plot_colors,loq_thresh=LOQ_sgRNA)

Plot_sumDyn_sgRNA_Naso <- Plot_sumDyn_sgRNA_Naso +
  ylab("Median sgRNA (copies/mL)") +
  ggtitle("Nasopharynx") + 
  theme(legend.position.inside = c(0.65,0.80),
        plot.title = element_text(face = "bold",color="black",size=13,hjust = 0)) +
  coord_cartesian(ylim=c(0.8*min(Viral_Load_sumDyn_sgRNA_Naso$Median),1E9),xlim=c(-0.25,30),expand = 0)
# ggsave(Plot_sumDyn_sgRNA_Naso,filename = paste(Figure_Folder,"Median_VL_sgRNA_Naso_dynamics_PAPER.png",sep="/"),dpi=300, height = 3,width=6)


## > Merged plots ####
Merged_Figure <- plot_grid(Plot_sumDyn_gRNA_Trachea,Plot_sumDyn_gRNA_Naso,Plot_sumDyn_sgRNA_Trachea,Plot_sumDyn_sgRNA_Naso,
                           labels=c("A","B","C","D"),
                           label_fontface = "bold",
                           ncol=2)
# ggsave(Merged_Figure,filename = paste(Figure_Folder,"Median_VL_dynamics_Merged_Figures_PAPER.png",sep="/"),height = 7,width=12, dpi=300)
# ---------------- #




# .................................... ----
# --- Calculation of VL dynamics characteristics (by group) ----
# .................................... ----
# For statistical tests performed on VL descriptors:
#  - For peak VL and AUC, normal distribution is considered (verified by Shapiro-Wilk normality test)
#        -> perform T-test for 2-2 mean comparisons (welch's or student t-test according to variance equality)
#        -> perform ANOVA for global mean comparisons
# - For time descriptors (time to peak, duration of clearance and acute stages), normality hypothesis is rejected 
#        -> perform Mann-whitney (= wilcoxon) test for 2-2 comparison 
#        -> perform Kruskal-Wallis for global comparison



Viral_Load_data_PostExpo <- subset(Viral_Load_data,SampleType %in% c("Nasal fluid","Tracheal fluid") & Time >=0)
Viral_Load_data_PostExpo <- Viral_Load_data_PostExpo[order(Viral_Load_data_PostExpo$Time),]

Overall_VL_descriptors <- data.frame(VLType=character(),SampleType=character(),
                                     Group=character(),AnimalID=character(),
                                     Descriptor=character(),Value=numeric())

## 1. Peak (value and time) of VL dynamics ----
### > Calculation of descriptors ----

Peak_VL <- ddply(.data=subset(Viral_Load_data_PostExpo,Censored == 0),.variables = .(VLType,SampleType,Group,SubjectID),summarize,
                 Peak = max(log10VLcens,na.rm=TRUE),
                 TimePeak = Time[which.max(log10VLcens)])

Overall_VL_descriptors <- rbind(Overall_VL_descriptors,melt(data=Peak_VL,id.vars = c("VLType","SampleType","Group","SubjectID"),variable.name = "Descriptor",value.name = "Value"))

Distribution_PeakVL <- ddply(.data=Peak_VL,.variables = .(VLType,SampleType,Group),summarize,
                             N = length(SubjectID),
                             Mean_Peak = mean(Peak,na.rm=TRUE),
                             ICMIN_Peak = quantile(Peak,na.rm=TRUE,probs=c(0.025)),
                             ICMAX_Peak = quantile(Peak,na.rm=TRUE,probs=c(0.975)),
                             
                             Mean_TimePeak = mean(TimePeak,na.rm=TRUE),
                             SD_TimePeak = sd(TimePeak,na.rm=TRUE),
                             ICMIN_TimePeak = quantile(TimePeak,na.rm=TRUE,probs=c(0.025)),
                             ICMAX_TimePeak = quantile(TimePeak,na.rm=TRUE,probs=c(0.975)))

### > Global Statistical test Comparison ----
Global_TestComparison_Peak <- expand.grid(VLType=c("gRNA","sgRNA"),SampleType=c("Tracheal fluid","Nasal fluid"),Descriptor=c("Peak","TimePeak"))
Global_TestComparison_Peak$Test <- ifelse(Global_TestComparison_Peak$Descriptor == "Peak","anova","ks")
Global_TestComparison_Peak$Pvalue <- sapply(1:nrow(Global_TestComparison_Peak),function(i) Global_TestComparison(data=subset(Peak_VL,VLType == Global_TestComparison_Peak$VLType[i] & SampleType == Global_TestComparison_Peak$SampleType[i]),
                                                                                                                 formula = formula(ifelse(Global_TestComparison_Peak$Descriptor[i] == "Peak",'Peak ~ Group','TimePeak~Group')),
                                                                                                                 test_name = ifelse(Global_TestComparison_Peak$Descriptor[i] == "Peak","anova","ks")))
Global_TestComparison_Peak$Adj.Pvalue <- p.adjust(Global_TestComparison_Peak$Pvalue,method="BH")


### > Two-sample comparison ----
TwoSample_test_Peak <- data.frame(gtools::combinations(n = 6, r = 2, repeats.allowed = F, v = as.character(unique(Peak_VL$Group))))
TwoSample_test_Peak <- rbind(data.frame(VLType = "gRNA",SampleType="Nasal fluid",Group1 = TwoSample_test_Peak$X1, Group2 = TwoSample_test_Peak$X2),
                             data.frame(VLType = "gRNA",SampleType="Tracheal fluid",Group1 = TwoSample_test_Peak$X1, Group2 = TwoSample_test_Peak$X2),
                             data.frame(VLType = "sgRNA",SampleType="Nasal fluid",Group1 = TwoSample_test_Peak$X1, Group2 = TwoSample_test_Peak$X2),
                             data.frame(VLType = "sgRNA",SampleType="Tracheal fluid",Group1 = TwoSample_test_Peak$X1, Group2 = TwoSample_test_Peak$X2))

TwoSample_test_Peak$Pvalue <- NA
TwoSample_test_TimePeak <- TwoSample_test_Peak

for(i in 1:nrow(TwoSample_test_Peak)){
  # i <- 1
  
  tmp_data_G1 <- subset(Peak_VL,VLType == TwoSample_test_Peak$VLType[i] & SampleType == TwoSample_test_Peak$SampleType[i] & Group == TwoSample_test_Peak$Group1[i] & !is.na(Peak))
  tmp_data_G2 <- subset(Peak_VL,VLType == TwoSample_test_Peak$VLType[i] & SampleType == TwoSample_test_Peak$SampleType[i] & Group == TwoSample_test_Peak$Group2[i] & !is.na(Peak))
  
  
  if(nrow(tmp_data_G1) > 1 & nrow(tmp_data_G2) > 1){
    TwoSample_test_Peak$Pvalue[i] <- TwoGroups_TestComparison(x=tmp_data_G1$Peak,y=tmp_data_G2$Peak,test_name = "ttest")
    TwoSample_test_TimePeak$Pvalue[i] <- TwoGroups_TestComparison(x=tmp_data_G1$TimePeak,y=tmp_data_G2$TimePeak,test_name="wilcoxon")
  }
}
TwoSample_test_Peak$Adj.Pvalue <- NA ; 
TwoSample_test_Peak$Adj.Pvalue[which(TwoSample_test_Peak$VLType == "gRNA" & TwoSample_test_Peak$SampleType == "Nasal fluid")] <- p.adjust(TwoSample_test_Peak$Pvalue[which(TwoSample_test_Peak$VLType == "gRNA" & TwoSample_test_Peak$SampleType == "Nasal fluid")],method="BH")
TwoSample_test_Peak$Adj.Pvalue[which(TwoSample_test_Peak$VLType == "gRNA" & TwoSample_test_Peak$SampleType == "Tracheal fluid")] <- p.adjust(TwoSample_test_Peak$Pvalue[which(TwoSample_test_Peak$VLType == "gRNA" & TwoSample_test_Peak$SampleType == "Tracheal fluid")],method="BH")
TwoSample_test_Peak$Adj.Pvalue[which(TwoSample_test_Peak$VLType == "sgRNA" & TwoSample_test_Peak$SampleType == "Nasal fluid")] <- p.adjust(TwoSample_test_Peak$Pvalue[which(TwoSample_test_Peak$VLType == "sgRNA" & TwoSample_test_Peak$SampleType == "Nasal fluid")],method="BH")
TwoSample_test_Peak$Adj.Pvalue[which(TwoSample_test_Peak$VLType == "sgRNA" & TwoSample_test_Peak$SampleType == "Tracheal fluid")] <- p.adjust(TwoSample_test_Peak$Pvalue[which(TwoSample_test_Peak$VLType == "sgRNA" & TwoSample_test_Peak$SampleType == "Tracheal fluid")],method="BH")
TwoSample_test_Peak$Significant <- 1*(TwoSample_test_Peak$Adj.Pvalue <= 0.05)

TwoSample_test_TimePeak$Adj.Pvalue <- NA
TwoSample_test_TimePeak$Adj.Pvalue[which(TwoSample_test_TimePeak$VLType == "gRNA" & TwoSample_test_TimePeak$SampleType == "Nasal fluid")] <- p.adjust(TwoSample_test_TimePeak$Pvalue[which(TwoSample_test_TimePeak$VLType == "gRNA" & TwoSample_test_TimePeak$SampleType == "Nasal fluid")],method="BH")
TwoSample_test_TimePeak$Adj.Pvalue[which(TwoSample_test_TimePeak$VLType == "gRNA" & TwoSample_test_TimePeak$SampleType == "Tracheal fluid")] <- p.adjust(TwoSample_test_TimePeak$Pvalue[which(TwoSample_test_TimePeak$VLType == "gRNA" & TwoSample_test_TimePeak$SampleType == "Tracheal fluid")],method="BH")
TwoSample_test_TimePeak$Adj.Pvalue[which(TwoSample_test_TimePeak$VLType == "sgRNA" & TwoSample_test_TimePeak$SampleType == "Nasal fluid")] <- p.adjust(TwoSample_test_TimePeak$Pvalue[which(TwoSample_test_TimePeak$VLType == "sgRNA" & TwoSample_test_TimePeak$SampleType == "Nasal fluid")],method="BH")
TwoSample_test_TimePeak$Adj.Pvalue[which(TwoSample_test_TimePeak$VLType == "sgRNA" & TwoSample_test_TimePeak$SampleType == "Tracheal fluid")] <- p.adjust(TwoSample_test_TimePeak$Pvalue[which(TwoSample_test_TimePeak$VLType == "sgRNA" & TwoSample_test_TimePeak$SampleType == "Tracheal fluid")],method="BH")
TwoSample_test_TimePeak$Significant <- 1*(TwoSample_test_TimePeak$Adj.Pvalue <= 0.05)


Distribution_PeakVL$Label_Peak <- with(Distribution_PeakVL,
                                       paste(round(Mean_Peak,digits = 2), " [",round(ICMIN_Peak,digits=2), " ; ",
                                             round(ICMAX_Peak,digits=2),"]",sep = ""))
Distribution_PeakVL$Label_TimePeak <- with(Distribution_PeakVL,
                                           paste(round(Mean_TimePeak,digits = 1), " [",round(ICMIN_TimePeak,digits=1), " ; ",
                                                 round(ICMAX_TimePeak,digits=1),"]",sep = ""))



## 2. AUC of VL dynamics ----
# Remarks: 
#    1) Censored data are set at the LOD value (i.e. AUC calculated above LOD)
#    2) To avoid any bias, we need to ensure that all NHPs are followed during the same time (the first 23 days)

### > Calculation of descriptors ----
Individual_AUC <- ddply(.data=subset(Viral_Load_data_PostExpo,Time <= 23),.variables = .(VLType,SampleType,Group,SubjectID),summarize,
                        AUC = AUC(x = Time,y = log10VLcens))

Overall_VL_descriptors <- rbind(Overall_VL_descriptors,melt(data=Individual_AUC,id.vars = c("VLType","SampleType","Group","SubjectID"),variable.name = "Descriptor",value.name = "Value"))


Distribution_AUC <- ddply(.data = Individual_AUC,.variables = .(VLType,SampleType,Group),summarize,
                          Mean_AUC = mean(AUC,na.rm = TRUE),
                          ICMIN_AUC = quantile(AUC,na.rm=TRUE,probs = c(0.025)),
                          ICMAX_AUC = quantile(AUC,na.rm=TRUE,probs = c(0.975)))

Distribution_AUC$Label <- with(Distribution_AUC,
                               paste(round(Mean_AUC,digits = 2), " [",round(ICMIN_AUC,digits=2), " ; ",
                                     round(ICMAX_AUC,digits=2),"]",sep = ""))


### > Global Statistical test Comparison ----
Global_TestComparison_AUC <- expand.grid(VLType=c("gRNA","sgRNA"),SampleType=c("Tracheal fluid","Nasal fluid"),Descriptor=c("AUC"))
Global_TestComparison_AUC$Test <- "anova"
Global_TestComparison_AUC$Pvalue <- sapply(1:nrow(Global_TestComparison_AUC),function(i) Global_TestComparison(data=subset(Individual_AUC,VLType == Global_TestComparison_AUC$VLType[i] & SampleType == Global_TestComparison_AUC$SampleType[i]),
                                                                                                               formula = formula(AUC ~ Group),
                                                                                                               test_name = Global_TestComparison_AUC$Test[i]))
Global_TestComparison_AUC$Adj.Pvalue <- p.adjust(Global_TestComparison_AUC$Pvalue,method="BH")

### > Two-sample comparison ----
TwoSample_test_AUC <- data.frame(gtools::combinations(n = 6, r = 2, repeats.allowed = F, v = as.character(unique(Individual_AUC$Group))))
TwoSample_test_AUC <- rbind(data.frame(VLType = "gRNA",SampleType="Nasal fluid",Group1 = TwoSample_test_AUC$X1, Group2 = TwoSample_test_AUC$X2),
                            data.frame(VLType = "gRNA",SampleType="Tracheal fluid",Group1 = TwoSample_test_AUC$X1, Group2 = TwoSample_test_AUC$X2),
                            data.frame(VLType = "sgRNA",SampleType="Nasal fluid",Group1 = TwoSample_test_AUC$X1, Group2 = TwoSample_test_AUC$X2),
                            data.frame(VLType = "sgRNA",SampleType="Tracheal fluid",Group1 = TwoSample_test_AUC$X1, Group2 = TwoSample_test_AUC$X2))

TwoSample_test_AUC$Pvalue <- NA
for(i in 1:nrow(TwoSample_test_AUC)){
  # i <- 34
  tmp_data_G1 <- subset(Individual_AUC,VLType == TwoSample_test_AUC$VLType[i] & SampleType == TwoSample_test_AUC$SampleType[i] & Group == TwoSample_test_AUC$Group1[i] & !is.na(AUC))
  tmp_data_G2 <- subset(Individual_AUC,VLType == TwoSample_test_AUC$VLType[i] & SampleType == TwoSample_test_AUC$SampleType[i] & Group == TwoSample_test_AUC$Group2[i] & !is.na(AUC))
  
  
  if(nrow(tmp_data_G1) > 1 & nrow(tmp_data_G2) > 1){
    TwoSample_test_AUC$Pvalue[i] <- TwoGroups_TestComparison(x=tmp_data_G1$AUC,y=tmp_data_G2$AUC,test_name = "ttest")
  }
}
TwoSample_test_AUC$Adj.Pvalue <- NA ; 
TwoSample_test_AUC$Adj.Pvalue[which(TwoSample_test_AUC$VLType == "gRNA" & TwoSample_test_AUC$SampleType == "Nasal fluid")] <- p.adjust(TwoSample_test_AUC$Pvalue[which(TwoSample_test_AUC$VLType == "gRNA" & TwoSample_test_AUC$SampleType == "Nasal fluid")],method="BH")
TwoSample_test_AUC$Adj.Pvalue[which(TwoSample_test_AUC$VLType == "gRNA" & TwoSample_test_AUC$SampleType == "Tracheal fluid")] <- p.adjust(TwoSample_test_AUC$Pvalue[which(TwoSample_test_AUC$VLType == "gRNA" & TwoSample_test_AUC$SampleType == "Tracheal fluid")],method="BH")
TwoSample_test_AUC$Adj.Pvalue[which(TwoSample_test_AUC$VLType == "sgRNA" & TwoSample_test_AUC$SampleType == "Nasal fluid")] <- p.adjust(TwoSample_test_AUC$Pvalue[which(TwoSample_test_AUC$VLType == "sgRNA" & TwoSample_test_AUC$SampleType == "Nasal fluid")],method="BH")
TwoSample_test_AUC$Adj.Pvalue[which(TwoSample_test_AUC$VLType == "sgRNA" & TwoSample_test_AUC$SampleType == "Tracheal fluid")] <- p.adjust(TwoSample_test_AUC$Pvalue[which(TwoSample_test_AUC$VLType == "sgRNA" & TwoSample_test_AUC$SampleType == "Tracheal fluid")],method="BH")
TwoSample_test_AUC$Significant <- 1*(TwoSample_test_AUC$Adj.Pvalue <= 0.05)




## 3. Duration of the clearance and acute stages ----
### > Calculation of descriptors ----
Individual_Stage <- unique(Viral_Load_data_PostExpo[,c("VLType","SampleType","Group","SubjectID")])
Individual_Stage$DurationClearanceStage <- NA ; Individual_Stage$DurationAcuteStage <- NA


for(i in 1:nrow(Individual_Stage)){
  
  # i <- 68
  # print(i)
  tmp_data <- subset(Viral_Load_data_PostExpo,VLType == Individual_Stage$VLType[i] & 
                       SampleType == Individual_Stage$SampleType[i] & SubjectID == Individual_Stage$SubjectID[i])
  tmp_peak <- subset(Peak_VL,VLType == Individual_Stage$VLType[i] & 
                       SampleType == Individual_Stage$SampleType[i] & SubjectID == Individual_Stage$SubjectID[i])
  
  if(nrow(tmp_peak) !=0){
    # Extraction of data related to clearance stage
    data_clearance_stage <- subset(tmp_data,Time > tmp_peak$TimePeak)
    firstUndetectableTime <- min(data_clearance_stage$Time[which(data_clearance_stage$Censored == 1)])
    Individual_Stage$DurationClearanceStage[i] <- firstUndetectableTime - tmp_peak$TimePeak
    
    # Extraction of data related to acute stage (detectable data)
    data_acute_stage <- subset(tmp_data,Censored == 0)
    Individual_Stage$DurationAcuteStage[i] <- max(data_acute_stage$Time) - min(data_acute_stage$Time)
  }
}

Overall_VL_descriptors <- rbind(Overall_VL_descriptors,melt(data=Individual_Stage,id.vars = c("VLType","SampleType","Group","SubjectID"),variable.name = "Descriptor",value.name = "Value"))


Distribution_Stages <- ddply(.data=Individual_Stage,.variables = .(VLType,SampleType,Group),summarize,
                             Mean_ClearanceStage = mean(DurationClearanceStage,na.rm=TRUE),
                             ICMIN_ClearanceStage = quantile(DurationClearanceStage,na.rm=TRUE,probs=c(0.025)),
                             ICMAX_ClearanceStage = quantile(DurationClearanceStage,na.rm=TRUE,probs=c(0.975)),
                             
                             Mean_AcuteStage = mean(DurationAcuteStage,na.rm=TRUE),
                             ICMIN_AcuteStage = quantile(DurationAcuteStage,na.rm=TRUE,probs=c(0.025)),
                             ICMAX_AcuteStage = quantile(DurationAcuteStage,na.rm=TRUE,probs=c(0.975)))

Distribution_Stages$Label_ClearanceStage <- with(Distribution_Stages,
                                                 paste(round(Mean_ClearanceStage,digits = 1), " [",round(ICMIN_ClearanceStage,digits=1), " ; ",
                                                       round(ICMAX_ClearanceStage,digits=1),"]",sep = ""))
Distribution_Stages$Label_AcuteStage <- with(Distribution_Stages,
                                             paste(round(Mean_AcuteStage,digits = 1), " [",round(ICMIN_AcuteStage,digits=1), " ; ",
                                                   round(ICMAX_AcuteStage,digits=1),"]",sep = ""))

### > Global Statistical test Comparison ----
Global_TestComparison_Stages <- expand.grid(VLType=c("gRNA","sgRNA"),SampleType=c("Tracheal fluid","Nasal fluid"),Descriptor=c("Clearance","Acute"))
Global_TestComparison_Stages$Test <- "ks"
Global_TestComparison_Stages$Pvalue <- sapply(1:nrow(Global_TestComparison_Stages),function(i) Global_TestComparison(data=subset(Individual_Stage,VLType == Global_TestComparison_Stages$VLType[i] & SampleType == Global_TestComparison_Stages$SampleType[i]),
                                                                                                                     formula = formula(ifelse(Global_TestComparison_Stages$Descriptor[i] == "Clearance",'DurationClearanceStage ~ Group','DurationAcuteStage~Group')),
                                                                                                                     test_name = Global_TestComparison_Stages$Test[i]))
Global_TestComparison_Stages$Adj.Pvalue <- p.adjust(Global_TestComparison_Stages$Pvalue,method="BH")



### > Two-sample comparison ----
TwoSample_test_ClearanceStage <- data.frame(gtools::combinations(n = 6, r = 2, repeats.allowed = F, v = as.character(unique(Individual_Stage$Group))))
TwoSample_test_ClearanceStage <- rbind(data.frame(VLType = "gRNA",SampleType="Nasal fluid",Group1 = TwoSample_test_ClearanceStage$X1, Group2 = TwoSample_test_ClearanceStage$X2),
                                       data.frame(VLType = "gRNA",SampleType="Tracheal fluid",Group1 = TwoSample_test_ClearanceStage$X1, Group2 = TwoSample_test_ClearanceStage$X2),
                                       data.frame(VLType = "sgRNA",SampleType="Nasal fluid",Group1 = TwoSample_test_ClearanceStage$X1, Group2 = TwoSample_test_ClearanceStage$X2),
                                       data.frame(VLType = "sgRNA",SampleType="Tracheal fluid",Group1 = TwoSample_test_ClearanceStage$X1, Group2 = TwoSample_test_ClearanceStage$X2))

TwoSample_test_ClearanceStage$Pvalue <- NA
TwoSample_test_AcuteStage <- TwoSample_test_ClearanceStage

for(i in 1:nrow(TwoSample_test_ClearanceStage)){
  # i <- 2
  
  tmp_data_G1 <- subset(Individual_Stage,VLType == TwoSample_test_ClearanceStage$VLType[i] & SampleType == TwoSample_test_ClearanceStage$SampleType[i] & Group == TwoSample_test_ClearanceStage$Group1[i] & !is.na(DurationClearanceStage))
  tmp_data_G2 <- subset(Individual_Stage,VLType == TwoSample_test_ClearanceStage$VLType[i] & SampleType == TwoSample_test_ClearanceStage$SampleType[i] & Group == TwoSample_test_ClearanceStage$Group2[i] & !is.na(DurationClearanceStage))
  
  
  if(nrow(tmp_data_G1) > 1 & nrow(tmp_data_G2) > 1){
    TwoSample_test_ClearanceStage$Pvalue[i] <- TwoGroups_TestComparison(x=tmp_data_G1$DurationClearanceStage,y=tmp_data_G2$DurationClearanceStage,test_name = "wilcoxon")
    TwoSample_test_AcuteStage$Pvalue[i] <- TwoGroups_TestComparison(x=tmp_data_G1$DurationAcuteStage,y=tmp_data_G2$DurationAcuteStage,test_name = "wilcoxon")
  }
}
TwoSample_test_ClearanceStage$Adj.Pvalue <- NA ; 
TwoSample_test_ClearanceStage$Adj.Pvalue[which(TwoSample_test_ClearanceStage$VLType == "gRNA" & TwoSample_test_ClearanceStage$SampleType == "Nasal fluid")] <- p.adjust(TwoSample_test_ClearanceStage$Pvalue[which(TwoSample_test_ClearanceStage$VLType == "gRNA" & TwoSample_test_ClearanceStage$SampleType == "Nasal fluid")],method="BH")
TwoSample_test_ClearanceStage$Adj.Pvalue[which(TwoSample_test_ClearanceStage$VLType == "gRNA" & TwoSample_test_ClearanceStage$SampleType == "Tracheal fluid")] <- p.adjust(TwoSample_test_ClearanceStage$Pvalue[which(TwoSample_test_ClearanceStage$VLType == "gRNA" & TwoSample_test_ClearanceStage$SampleType == "Tracheal fluid")],method="BH")
TwoSample_test_ClearanceStage$Adj.Pvalue[which(TwoSample_test_ClearanceStage$VLType == "sgRNA" & TwoSample_test_ClearanceStage$SampleType == "Nasal fluid")] <- p.adjust(TwoSample_test_ClearanceStage$Pvalue[which(TwoSample_test_ClearanceStage$VLType == "sgRNA" & TwoSample_test_ClearanceStage$SampleType == "Nasal fluid")],method="BH")
TwoSample_test_ClearanceStage$Adj.Pvalue[which(TwoSample_test_ClearanceStage$VLType == "sgRNA" & TwoSample_test_ClearanceStage$SampleType == "Tracheal fluid")] <- p.adjust(TwoSample_test_ClearanceStage$Pvalue[which(TwoSample_test_ClearanceStage$VLType == "sgRNA" & TwoSample_test_ClearanceStage$SampleType == "Tracheal fluid")],method="BH")
TwoSample_test_ClearanceStage$Significant <- 1*(TwoSample_test_ClearanceStage$Adj.Pvalue <= 0.05)

TwoSample_test_AcuteStage$Adj.Pvalue <- NA
TwoSample_test_AcuteStage$Adj.Pvalue[which(TwoSample_test_AcuteStage$VLType == "gRNA" & TwoSample_test_AcuteStage$SampleType == "Nasal fluid")] <- p.adjust(TwoSample_test_AcuteStage$Pvalue[which(TwoSample_test_AcuteStage$VLType == "gRNA" & TwoSample_test_AcuteStage$SampleType == "Nasal fluid")],method="BH")
TwoSample_test_AcuteStage$Adj.Pvalue[which(TwoSample_test_AcuteStage$VLType == "gRNA" & TwoSample_test_AcuteStage$SampleType == "Tracheal fluid")] <- p.adjust(TwoSample_test_AcuteStage$Pvalue[which(TwoSample_test_AcuteStage$VLType == "gRNA" & TwoSample_test_AcuteStage$SampleType == "Tracheal fluid")],method="BH")
TwoSample_test_AcuteStage$Adj.Pvalue[which(TwoSample_test_AcuteStage$VLType == "sgRNA" & TwoSample_test_AcuteStage$SampleType == "Nasal fluid")] <- p.adjust(TwoSample_test_AcuteStage$Pvalue[which(TwoSample_test_AcuteStage$VLType == "sgRNA" & TwoSample_test_AcuteStage$SampleType == "Nasal fluid")],method="BH")
TwoSample_test_AcuteStage$Adj.Pvalue[which(TwoSample_test_AcuteStage$VLType == "sgRNA" & TwoSample_test_AcuteStage$SampleType == "Tracheal fluid")] <- p.adjust(TwoSample_test_AcuteStage$Pvalue[which(TwoSample_test_AcuteStage$VLType == "sgRNA" & TwoSample_test_AcuteStage$SampleType == "Tracheal fluid")],method="BH")
TwoSample_test_AcuteStage$Significant <- 1*(TwoSample_test_AcuteStage$Adj.Pvalue <= 0.05)







## 4. Plot of results ----
Plot_colors <- Plot_colors_Groups()

Overall_VL_descriptors$Group <- factor(Overall_VL_descriptors$Group,levels = names(Plot_colors))
Overall_VL_descriptors$Descriptor <- as.character(Overall_VL_descriptors$Descriptor)
Overall_VL_descriptors$Descriptor[which(Overall_VL_descriptors$Descriptor == "DurationClearanceStage")] <- "Clearance"
Overall_VL_descriptors$Descriptor[which(Overall_VL_descriptors$Descriptor == "DurationAcuteStage")] <- "Acute"
Overall_VL_descriptors$Descriptor <- factor(Overall_VL_descriptors$Descriptor,levels=c("TimePeak","Peak","AUC","Clearance","Acute"))

# Calculation of the distribution of each descriptor
Distribution_Descriptors <- ddply(.data = Overall_VL_descriptors,.variables = .(VLType,SampleType,Group,Descriptor),summarize,
                                  Mean = mean(Value,na.rm=TRUE),
                                  Median = median(Value,na.rm=TRUE),
                                  ICMIN=quantile(Value,0.025,na.rm=TRUE),
                                  IC25=quantile(Value,0.25,na.rm=TRUE),
                                  IC75=quantile(Value,0.75,na.rm=TRUE),
                                  ICMAX=quantile(Value,0.975,na.rm=TRUE))


Global_TestComparison <- rbind(Global_TestComparison_Peak,Global_TestComparison_AUC,Global_TestComparison_Stages)
Global_TestComparison$Significant <- with(Global_TestComparison,1*(Adj.Pvalue<=0.05))
Global_TestComparison$Label <- with(Global_TestComparison,
                                    case_when(Test == "anova" & Adj.Pvalue < 1e-3 ~ "ANOVA pvalue < 0.001",
                                              Test == "ks" & Adj.Pvalue < 1e-3 ~ "KS pvalue < 0.001",
                                              Test == "anova" & Adj.Pvalue >= 1e-3 ~ paste("ANOVA pvalue = ",round(Adj.Pvalue,digits=3),sep=""),
                                              Test == "ks" & Adj.Pvalue >= 1e-3 ~ paste("KS pvalue = ",round(Adj.Pvalue,digits=3),sep="")))
Global_TestComparison$x <- with(Global_TestComparison,
                                case_when(Descriptor == "TimePeak" & VLType == "gRNA" ~ 1,
                                          Descriptor == "TimePeak" & VLType == "sgRNA" ~ 1,
                                          Descriptor == "Peak" & VLType == "gRNA" ~ 3,
                                          Descriptor == "Peak" & VLType == "sgRNA" ~ 3,
                                          Descriptor == "AUC" & VLType == "gRNA" ~ 55,
                                          Descriptor == "AUC" & VLType == "sgRNA" ~ 55,
                                          Descriptor %in% c("Clearance","Acute") ~ 0,
                                          TRUE ~ 0))
Global_TestComparison$y <- 7.2

Xlimits_Labels <- unique(Global_TestComparison[,c("VLType","SampleType","Descriptor")])
Xlimits_Labels$xmin <- with(Xlimits_Labels,
                            case_when(Descriptor == "TimePeak" ~ 1,
                                      Descriptor == "Peak" ~ 3,
                                      Descriptor == "AUC"  ~ 55,
                                      Descriptor %in% c("Clearance","Acute") ~ 0))
Xlimits_Labels$xmax <- with(Xlimits_Labels,
                            case_when(Descriptor == "TimePeak" & SampleType == "Tracheal fluid" ~ 3,
                                      Descriptor == "TimePeak" & SampleType == "Nasal fluid" ~ 7,
                                      Descriptor == "Peak" ~ 9,
                                      Descriptor == "AUC" & SampleType == "Tracheal fluid" ~ 110,
                                      Descriptor == "AUC" & SampleType == "Nasal fluid" ~ 130,
                                      Descriptor == "Clearance" & SampleType == "Tracheal fluid" ~ 20,
                                      Descriptor == "Clearance" & SampleType == "Nasal fluid" ~ 30,
                                      Descriptor == "Acute" & SampleType == "Tracheal fluid" ~ 13,
                                      Descriptor == "Acute" & SampleType == "Nasal fluid" ~ 30))
Xlimits_Labels$y <- 7.2

Labels_Descriptors <- c("Time to peak viral load \n [in days]",
                        "Peak viral load \n [in log10 copies/mL]",
                        "Area under viral load curve \n [in log10 copies.day/mL]",
                        "Duration clearance stage \n [in days]",
                        "Duration acute stage \n [in days]")
names(Labels_Descriptors) <- c("TimePeak","Peak","AUC","Clearance","Acute")



### > gRNA and sgRNA in Trachea ----
GeomPointRange_gRNA_Trachea <- ggplot(data=subset(Distribution_Descriptors,VLType == "gRNA" & SampleType == "Tracheal fluid")) + 
  geom_pointrange(aes(x=Mean,xmin=ICMIN,xmax=ICMAX,y=Group,color=Group),shape=21,stroke=1.5,cex=0.60,fill="white",linewidth=1.0,show.legend = FALSE) + 
  facet_rep_grid(.~Descriptor,scales="free",repeat.tick.labels=FALSE,labeller=labeller(Descriptor = Labels_Descriptors)) + 
  geom_text(data=subset(Global_TestComparison,VLType == "gRNA" & SampleType == "Tracheal fluid"),
            aes(x=x,y=y,label=Label),hjust = 0,size = 3.5,fontface = "plain") + 
  geom_text(data=Xlimits_Labels,aes(x=xmin,y=y,label = "."),color="white") +
  geom_text(data=Xlimits_Labels,aes(x=xmax,y=y,label = "."),color="white") +
  
  scale_color_manual(values = Plot_colors) + 
  coord_cartesian(ylim=c(1,7)) +
  ggtitle("Genomic RNA in Trachea") + 
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


GeomPointRange_sgRNA_Trachea <- ggplot(data=subset(Distribution_Descriptors,VLType == "sgRNA" & SampleType == "Tracheal fluid")) + 
  geom_pointrange(aes(x=Mean,xmin=ICMIN,xmax=ICMAX,y=Group,color=Group),shape=21,stroke=1.5,cex=0.60,fill="white",linewidth=1.0,show.legend = FALSE) + 
  facet_rep_grid(.~Descriptor,scales="free",repeat.tick.labels=FALSE,labeller=labeller(Descriptor = Labels_Descriptors)) + 
  geom_text(data=subset(Global_TestComparison,VLType == "sgRNA" & SampleType == "Tracheal fluid"),
            aes(x=x,y=y,label=Label),hjust = 0,size = 3.5,fontface = "plain") + 
  geom_text(data=Xlimits_Labels,aes(x=xmin,y=y,label = "."),color="white") +
  geom_text(data=Xlimits_Labels,aes(x=xmax,y=y,label = "."),color="white") +
  
  scale_color_manual(values = Plot_colors) + 
  coord_cartesian(ylim=c(1,7)) +
  ggtitle("Subgenomic RNA in Trachea") + 
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

# Merge of gRNA and sgRNA plots
GeomPointRange_VL_Trachea <- cowplot::plot_grid(GeomPointRange_gRNA_Trachea,GeomPointRange_sgRNA_Trachea,nrow = 2,align = "h",labels = "AUTO")
# ggsave(plot=GeomPointRange_VL_Trachea,filename = paste(Figure_Folder,"VLtrachea_Descriptors_PAPER.png",sep="/"),height = 6,width=12,dpi=300)



### > gRNA and sgRNA in Nasopharynx ----
GeomPointRange_gRNA_Naso <- ggplot(data=subset(Distribution_Descriptors,VLType == "gRNA" & SampleType == "Nasal fluid")) + 
  geom_pointrange(aes(x=Mean,xmin=ICMIN,xmax=ICMAX,y=Group,color=Group),shape=21,stroke=1.5,cex=0.60,fill="white",linewidth=1.0,show.legend = FALSE) + 
  facet_rep_grid(.~Descriptor,scales="free",repeat.tick.labels=FALSE,labeller=labeller(Descriptor = Labels_Descriptors)) + 
  geom_text(data=subset(Global_TestComparison,VLType == "gRNA" & SampleType == "Nasal fluid"),
            aes(x=x,y=y,label=Label),hjust = 0,size = 3.5,fontface = "plain") + 
  geom_text(data=Xlimits_Labels,aes(x=xmin,y=y,label = "."),color="white") +
  geom_text(data=Xlimits_Labels,aes(x=xmax,y=y,label = "."),color="white") +
  
  scale_color_manual(values = Plot_colors) + 
  coord_cartesian(ylim=c(1,7)) +
  ggtitle("Genomic RNA in Nasopharynx") + 
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


GeomPointRange_sgRNA_Naso <- ggplot(data=subset(Distribution_Descriptors,VLType == "sgRNA" & SampleType == "Nasal fluid")) + 
  geom_pointrange(aes(x=Mean,xmin=ICMIN,xmax=ICMAX,y=Group,color=Group),shape=21,stroke=1.5,cex=0.60,fill="white",linewidth=1.0,show.legend = FALSE) + 
  facet_rep_grid(.~Descriptor,scales="free",repeat.tick.labels=FALSE,labeller=labeller(Descriptor = Labels_Descriptors)) + 
  geom_text(data=subset(Global_TestComparison,VLType == "sgRNA" & SampleType == "Nasal fluid"),
            aes(x=x,y=y,label=Label),hjust = 0,size = 3.5,fontface = "plain") + 
  geom_text(data=Xlimits_Labels,aes(x=xmin,y=y,label = "."),color="white") +
  geom_text(data=Xlimits_Labels,aes(x=xmax,y=y,label = "."),color="white") +
  
  scale_color_manual(values = Plot_colors) + 
  coord_cartesian(ylim=c(1,7)) +
  ggtitle("Subgenomic RNA in Nasopharynx") + 
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

# Merge of gRNA and sgRNA plots
GeomPointRange_VL_Naso <- cowplot::plot_grid(GeomPointRange_gRNA_Naso,GeomPointRange_sgRNA_Naso,nrow = 2,align = "h",labels = "AUTO")
# ggsave(plot=GeomPointRange_VL_Naso,filename = paste(Figure_Folder,"VLnaso_Descriptors_PAPER.png",sep="/"),height = 6,width=12,dpi=300)
# ----------- #






