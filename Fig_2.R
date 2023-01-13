library(tidyverse)
library(lubridate)
library(cowplot)

# Purpose: Plot ELISA ODs using box and jitter plot in Figure 2 by vaccination and infection history
# Note 1: Vaccination history data cannot be shared pursuant to confidentiality Undertaking to the Department of Health, Government of the HKSAR. 
# Note 2: Serology data are available from the authors upon reasonable request.

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Read data
confirmed_vax_data <- read_csv("confirmed_vax_data.csv")

universal_thresholds_data <- read_csv ("lab_data_universal_thresholds.csv")
NCTD_OD_limit <- as.numeric(universal_thresholds_data%>%filter(antigen=="CTD") %>% select(`threshold`))
ORF8_OD_limit <- as.numeric(universal_thresholds_data%>%filter(antigen=="ORF8") %>% select(`threshold`))

confirmed_vax_data <- confirmed_vax_data %>% filter(vaccine_category=="BNT"|vaccine_category=="Sinovac"|vaccine_category=="never") %>% filter (!is.na(age_group))
confirmed_vax_data$vaccine_category <- as.factor(confirmed_vax_data$vaccine_category)
confirmed_vax_data$no_of_vaccines <- as.factor(confirmed_vax_data$no_of_vaccines)
confirmed_vax_data$previous_confirmed_infection_self <- as.factor(confirmed_vax_data$previous_confirmed_infection_self)

confirmed_vax_data$age_group <- as.factor(confirmed_vax_data$age_group)

levels(confirmed_vax_data$previous_confirmed_infection_self) <- list("No or unknown"=c("Choose not to answer","N",NA),
                                                                     "Yes"=c("Y"))

confirmed_vax_data$previous_confirmed_infection_self[is.na(confirmed_vax_data$previous_confirmed_infection_self)] <- "No or unknown"

# Plot unvaccinated ELISA ODs by age group and self-reported infection history
no_vax <- confirmed_vax_data %>% filter(no_of_vaccines=="0")
no_vax_plot <- no_vax %>% ggplot(aes(x=age_group,y=CTD_protein_OD))+
  geom_jitter(colour="#8491B499",alpha=0.3)+
  geom_boxplot(outlier.shape=NA,alpha=0)+
  geom_hline(yintercept = NCTD_OD_limit,colour="#3C5488FF",linetype="dashed",size=1)+
  labs(x="Age group",y="N-CTD OD")+
  facet_wrap(.~previous_confirmed_infection_self)+
  theme_classic()+
  theme(strip.text.x = element_text(size = 12))


no_vax_plot

# Plot homologous BNT-vaccinated ELISA ODs by age group and self-reported infection history
BNT_vax <- confirmed_vax_data %>% filter(vaccine_category=="BNT"&no_of_vaccines!="0")
BNT_vax_plot <- BNT_vax %>% ggplot(aes(x=age_group,y=CTD_protein_OD))+
  geom_jitter(colour="#8491B499",alpha=0.3)+
  geom_boxplot(outlier.shape=NA,alpha=0)+
  geom_hline(yintercept = NCTD_OD_limit,colour="#3C5488FF",linetype="dashed",size=1)+
  labs(x="Age group",y="N-CTD OD")+
  facet_wrap(.~previous_confirmed_infection_self)+
  theme_classic()+
  theme(strip.text.x = element_text(size = 12))

BNT_vax_plot

# Plot homologous CoronaVac-vaccinated ELISA ODs by age group and self-reported infection history
C_vax <- confirmed_vax_data %>% filter(vaccine_category=="Sinovac"&no_of_vaccines!="0")
C_vax_plot <- C_vax %>% ggplot(aes(x=age_group,y=ORF8_OD))+
  geom_jitter(colour="#8491B499",alpha=0.3)+
  geom_boxplot(outlier.shape=NA,alpha=0)+
  geom_hline(yintercept = NCTD_OD_limit,colour="#3C5488FF",linetype="dashed",size=1)+
  labs(x="Age group",y="ORF8 OD")+
  facet_wrap(.~previous_confirmed_infection_self)+
  theme_classic()+
  theme(strip.text.x = element_text(size = 12))


C_vax_plot

all_vax_plot <- plot_grid(no_vax_plot+theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_text(size=10),axis.text.x=element_text(size=10),axis.text.y=element_text(size=10)),
                          BNT_vax_plot+theme(legend.position="none",axis.title.x = element_blank(),axis.title.y = element_text(size=10),axis.text.x=element_text(size=10),axis.text.y=element_text(size=10)),
                          C_vax_plot+theme(legend.position="none",axis.title.x=element_text(size=10),axis.title.y = element_text(size=10),axis.text.x=element_text(size=10),axis.text.y=element_text(size=10)),
                          nrow=3,
                          rel_heights=c(1,1,1),
                          labels="AUTO",
                          label_x=0,
                          label_y=1)

ggsave(file="Extended_data_fig_4.tiff",all_vax_plot,width=8,height=8,dpi=1000,compression="lzw")
all_vax_plot
dev.off()
