library(pROC)
library(tidyverse)
library(cowplot)
library(xlsx)

# Purpose: To derive bootstrapped ROC curves using package pROC (Robin et al, BMC Bioinformatics 2011) 

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

threshold = c(seq(0,1,by=0.001))

# Read negative and positive control data

neg_controls_unvax <- read_csv("CTD_neg_controls_unvax.csv") %>% select(CTD.OD,D)
neg_controls_BNT <- read_csv("CTD_neg_controls_BNT.csv") %>% select(CTD.OD,D)
neg_controls_C <- read_csv("ORF8_neg_controls_C.csv") %>% select(ORF8.OD,D)

pos_controls_unvax <- read_csv("CTD_pos_controls_unvax.csv")
pos_controls_BNT <- read_csv("CTD_pos_controls_BNT.csv")
pos_controls_C <- read_csv("ORF8_pos_controls_C.csv")

performance <- read_csv("assay_performance.csv",col_names=c("Median","LCI","UCI"))
performance[c(4:6),] <- 1-performance[c(4:6),]
performance$test <- c("Sens_unvax","Sens_BNT","Sens_C","Spec_unvax","Spec_BNT","Spec_C")

controls_unvax <- rbind(neg_controls_unvax,pos_controls_unvax)
controls_BNT <- rbind(neg_controls_BNT,pos_controls_BNT)
controls_C <- rbind(neg_controls_C,pos_controls_C)

controls_BNT_unvax <- rbind(controls_unvax,controls_BNT)


# Calculate and plot ROC for unvaccinated individuals

roc_unvax.list <- roc(controls_unvax$D~controls_unvax$CTD.OD)
ci_unvax.list <- ci.se(roc_unvax.list,specificities = seq(0,1,l=100))
ci_unvax.df <- data.frame(x=1-as.numeric(rownames(ci_unvax.list)),ymin=ci_unvax.list[,1],ymax=ci_unvax.list[,3])
label_unvax <- paste(paste("Sens: ",round(performance$`Median`[performance$test=="Sens_unvax"],2)," (",round(performance$`LCI`[performance$test=="Sens_unvax"],2),",",round(performance$`UCI`[performance$test=="Sens_unvax"],2),")",sep=""),
                     paste("Spec: ",round(1-performance$`Median`[performance$test=="Spec_unvax"],2)," (",round(1-performance$`LCI`[performance$test=="Spec_unvax"],2),",",round(1-performance$`UCI`[performance$test=="Spec_unvax"],2),")",sep=""),
                     paste("AUC:",round(auc(roc_unvax.list),4)),sep="\n")
roc_unvax <- ggroc(roc_unvax.list,legacy.axes = TRUE)+theme_minimal()+geom_abline(slope=1,intercept=0,linetype="dashed",alpha=0.7,color="grey")+coord_equal()+
  geom_ribbon(data=ci_unvax.df,aes(x=x,ymin=ymin,ymax=ymax),inherit.aes = FALSE,fill="#56B4E9",alpha=0.5)+
  geom_point(aes(x=performance$`Median`[performance$test=="Spec_unvax"],y=performance$`Median`[performance$test=="Sens_unvax"]),colour="#E64B35FF")+
  geom_errorbar(aes(x=performance$`Median`[performance$test=="Spec_unvax"],ymin=performance$`LCI`[performance$test=="Sens_unvax"],ymax=performance$`UCI`[performance$test=="Sens_unvax"]),position=position_dodge(0.9),width=0.01,colour="#E64B35FF")+
  geom_errorbarh(aes(y=performance$`Median`[performance$test=="Sens_unvax"],xmin=performance$`UCI`[performance$test=="Spec_unvax"],xmax=performance$`LCI`[performance$test=="Spec_unvax"]),position=position_dodge(0.9),height=0.01,colour="#E64B35FF")+
  xlab("% false positive (1-specificity)")+
  ylab("% true positive (sensitivity)")+
  ggtitle(paste0("N-CTD:Unvaccinated\n","Pos:",nrow(pos_controls_unvax),", Neg:",nrow(neg_controls_unvax)))+
  theme(plot.title=element_text(face="bold",size=10))+
  geom_label(x=0.65,y=0.2,label=label_unvax,size=3)


roc_unvax
thresholds_unvax <- ci.thresholds(roc_unvax.list,thresholds=threshold)
tpr_unvax <- data.table::transpose(coords(roc_unvax.list,seq(0,1,0.05)))[2,]


# Calculate and plot ROC for homologous BNT-vaccinated individuals
roc_BNT.list <- roc(controls_BNT$D~controls_BNT$CTD.OD)
ci_BNT.list <- ci.se(roc_BNT.list,specificities = seq(0,1,l=100))
ci_BNT.df <- data.frame(x=1-as.numeric(rownames(ci_BNT.list)),ymin=ci_BNT.list[,1],ymax=ci_BNT.list[,3])
label_BNT <- paste(paste("Sens: ",round(performance$`Median`[performance$test=="Sens_BNT"],2)," (",round(performance$`LCI`[performance$test=="Sens_BNT"],2),",",round(performance$`UCI`[performance$test=="Sens_BNT"],2),")",sep=""),
                   paste("Spec: ",round(1-performance$`Median`[performance$test=="Spec_BNT"],2)," (",round(1-performance$`LCI`[performance$test=="Spec_BNT"],2),",",round(1-performance$`UCI`[performance$test=="Spec_BNT"],2),")",sep=""),
                   paste("AUC:",round(auc(roc_BNT.list),4)),sep="\n")
roc_BNT <- ggroc(roc_BNT.list,legacy.axes = TRUE)+theme_minimal()+geom_abline(slope=1,intercept=0,linetype="dashed",alpha=0.7,color="grey")+coord_equal()+
  geom_ribbon(data=ci_BNT.df,aes(x=x,ymin=ymin,ymax=ymax),inherit.aes = FALSE,fill="#56B4E9",alpha=0.5)+
  geom_point(aes(x=performance$`Median`[performance$test=="Spec_BNT"],y=performance$`Median`[performance$test=="Sens_BNT"]),colour="#E64B35FF")+
  geom_errorbar(aes(x=performance$`Median`[performance$test=="Spec_BNT"],ymin=performance$`LCI`[performance$test=="Sens_BNT"],ymax=performance$`UCI`[performance$test=="Sens_BNT"]),position=position_dodge(0.9),width=0.01,colour="#E64B35FF")+
  geom_errorbarh(aes(y=performance$`Median`[performance$test=="Sens_BNT"],xmin=performance$`UCI`[performance$test=="Spec_BNT"],xmax=performance$`LCI`[performance$test=="Spec_BNT"]),position=position_dodge(0.9),height=0.01,colour="#E64B35FF")+
  xlab("% false positive (1-specificity)")+
  ylab("% true positive (sensitivity)")+
  ggtitle(paste0("N-CTD:BNT\n","Pos:",nrow(pos_controls_BNT),", Neg:",nrow(neg_controls_BNT)))+
  theme(plot.title=element_text(face="bold",size=10))+
  geom_label(x=0.65,y=0.2,label=label_BNT,size=3)

roc_BNT
thresholds_BNT <- ci.thresholds(roc_BNT.list,thresholds=threshold)
tpr_BNT <- data.table::transpose(coords(roc_BNT.list,seq(0,1,0.05)))[2,]


# Calculate and plot ROC for homologous CoronaVac-vaccinated individuals
roc_C.list <- roc(controls_C$D~controls_C$ORF8.OD)
ci_C.list <- ci.se(roc_C.list,specificities = seq(0,1,l=100))
ci_C.df <- data.frame(x=1-as.numeric(rownames(ci_C.list)),ymin=ci_C.list[,1],ymax=ci_C.list[,3])
label_C <- paste(paste("Sens: ",round(performance$`Median`[performance$test=="Sens_C"],2)," (",round(performance$`LCI`[performance$test=="Sens_C"],2),",",round(performance$`UCI`[performance$test=="Sens_C"],2),")",sep=""),
                 paste("Spec: ",round(1-performance$`Median`[performance$test=="Spec_C"],2)," (",round(1-performance$`LCI`[performance$test=="Spec_C"],2),",",round(1-performance$`UCI`[performance$test=="Spec_C"],2),")",sep=""),
                 paste("AUC:",round(auc(roc_C.list),4)),sep="\n")
roc_C <- ggroc(roc_C.list,legacy.axes = TRUE)+theme_minimal()+geom_abline(slope=1,intercept=0,linetype="dashed",alpha=0.7,color="grey")+coord_equal()+
  geom_ribbon(data=ci_C.df,aes(x=x,ymin=ymin,ymax=ymax),inherit.aes = FALSE,fill="#56B4E9",alpha=0.5)+
  geom_point(aes(x=performance$`Median`[performance$test=="Spec_C"],y=performance$`Median`[performance$test=="Sens_C"]),colour="#E64B35FF")+
  geom_errorbar(aes(x=performance$`Median`[performance$test=="Spec_C"],ymin=performance$`LCI`[performance$test=="Sens_C"],ymax=performance$`UCI`[performance$test=="Sens_C"]),position=position_dodge(0.9),width=0.01,colour="#E64B35FF")+
  geom_errorbarh(aes(y=performance$`Median`[performance$test=="Sens_C"],xmin=performance$`UCI`[performance$test=="Spec_C"],xmax=performance$`LCI`[performance$test=="Spec_C"]),position=position_dodge(0.9),height=0.01,colour="#E64B35FF")+
  xlab("% false positive (1-specificity)")+
  ylab("% true positive (sensitivity)")+
  ggtitle(paste0("ORF8:CoronaVac\n","Pos:",nrow(pos_controls_C),", Neg:",nrow(neg_controls_C)))+
  theme(plot.title=element_text(face="bold",size=10))+
  geom_label(x=0.65,y=0.2,label=label_C,size=3)

roc_C
thresholds_C <- ci.thresholds(roc_C.list,thresholds=threshold)
tpr_C <- data.table::transpose(coords(roc_C.list,seq(0,1,0.05)))[2,]



roc_list <- list(
  unvax = roc_unvax,
  BNT = roc_BNT,
  sinovac = roc_C)

full_plot <- plot_grid(roc_unvax,roc_BNT,roc_C,ncol=3,nrow=1,labels="AUTO")

ggsave("Extended_Data_Fig_3.tiff",full_plot,width=9,height=3, dpi=1000, compression="lzw")

