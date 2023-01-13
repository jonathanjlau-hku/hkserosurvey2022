library (tidyverse)
library (forcats)
library(cowplot)
library(lubridate)

rm(list=ls())
# dev.off()

## Purpose: To create Extended Data Figure 4 of estimates of VE and waning using 
## only samples collected on or before June 15, 2022

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
startDate = ymd("2022-01-01")

# VE plots (Extended Data Figure 4)

VE_BB_prctile_7delay <- as.data.frame(cbind(rep("BNT",201),rep("2 doses",201),rep("7 days",201),read_csv("./VE/VE_BB_prctile_7delay_age.csv",col_names = FALSE)))
VE_BB_prctile_14delay <- as.data.frame(cbind(rep("BNT",201),rep("2 doses",201),rep("14 days",201),read_csv("./VE/VE_BB_prctile_14delay_age.csv",col_names = FALSE)))
VE_BB_prctile_21delay <- as.data.frame(cbind(rep("BNT",201),rep("2 doses",201),rep("21 days",201),read_csv("./VE/VE_BB_prctile_21delay_age.csv",col_names = FALSE)))

VE_BBB_prctile_7delay <- as.data.frame(cbind(rep("BNT",201),rep("3 doses",201),rep("7 days",201),read_csv("./VE/VE_BBB_prctile_7delay_age.csv",col_names = FALSE)))
VE_BBB_prctile_14delay <- as.data.frame(cbind(rep("BNT",201),rep("3 doses",201),rep("14 days",201),read_csv("./VE/VE_BBB_prctile_14delay_age.csv",col_names = FALSE)))
VE_BBB_prctile_21delay <- as.data.frame(cbind(rep("BNT",201),rep("3 doses",201),rep("21 days",201),read_csv("./VE/VE_BBB_prctile_21delay_age.csv",col_names = FALSE)))

VE_BBBB_prctile_7delay <- as.data.frame(cbind(rep("BNT",201),rep("4 doses",201),rep("7 days",201),read_csv("./VE/VE_BBBB_prctile_7delay_age.csv",col_names = FALSE)))
VE_BBBB_prctile_14delay <- as.data.frame(cbind(rep("BNT",201),rep("4 doses",201),rep("14 days",201),read_csv("./VE/VE_BBBB_prctile_14delay_age.csv",col_names = FALSE)))
VE_BBBB_prctile_21delay <- as.data.frame(cbind(rep("BNT",201),rep("4 doses",201),rep("21 days",201),read_csv("./VE/VE_BBBB_prctile_21delay_age.csv",col_names = FALSE)))

VE_CC_prctile_7delay <- as.data.frame(cbind(rep("CoronaVac",201),rep("2 doses",201),rep("7 days",201),read_csv("./VE/VE_CC_prctile_7delay_age.csv",col_names = FALSE)))
VE_CC_prctile_14delay <- as.data.frame(cbind(rep("CoronaVac",201),rep("2 doses",201),rep("14 days",201),read_csv("./VE/VE_CC_prctile_14delay_age.csv",col_names = FALSE)))
VE_CC_prctile_21delay <- as.data.frame(cbind(rep("CoronaVac",201),rep("2 doses",201),rep("21 days",201),read_csv("./VE/VE_CC_prctile_21delay_age.csv",col_names = FALSE)))

VE_CCC_prctile_7delay <- as.data.frame(cbind(rep("CoronaVac",201),rep("3 doses",201),rep("7 days",201),read_csv("./VE/VE_CCC_prctile_7delay_age.csv",col_names = FALSE)))
VE_CCC_prctile_14delay <- as.data.frame(cbind(rep("CoronaVac",201),rep("3 doses",201),rep("14 days",201),read_csv("./VE/VE_CCC_prctile_14delay_age.csv",col_names = FALSE)))
VE_CCC_prctile_21delay <- as.data.frame(cbind(rep("CoronaVac",201),rep("3 doses",201),rep("21 days",201),read_csv("./VE/VE_CCC_prctile_21delay_age.csv",col_names = FALSE)))

VE_CCCC_prctile_7delay <- as.data.frame(cbind(rep("CoronaVac",201),rep("4 doses",201),rep("7 days",201),read_csv("./VE/VE_CCCC_prctile_7delay_age.csv",col_names = FALSE)))
VE_CCCC_prctile_14delay <- as.data.frame(cbind(rep("CoronaVac",201),rep("4 doses",201),rep("14 days",201),read_csv("./VE/VE_CCCC_prctile_14delay_age.csv",col_names = FALSE)))
VE_CCCC_prctile_21delay <- as.data.frame(cbind(rep("CoronaVac",201),rep("4 doses",201),rep("21 days",201),read_csv("./VE/VE_CCCC_prctile_21delay_age.csv",col_names = FALSE)))


for (df in c("VE_BB_prctile_7delay","VE_BB_prctile_14delay","VE_BB_prctile_21delay","VE_BBB_prctile_7delay","VE_BBB_prctile_14delay","VE_BBB_prctile_21delay","VE_BBBB_prctile_7delay","VE_BBBB_prctile_14delay","VE_BBBB_prctile_21delay",
             "VE_CC_prctile_7delay","VE_CC_prctile_14delay","VE_CC_prctile_21delay","VE_CCC_prctile_7delay","VE_CCC_prctile_14delay","VE_CCC_prctile_21delay","VE_CCCC_prctile_7delay","VE_CCCC_prctile_14delay","VE_CCCC_prctile_21delay"))
  data.table::setnames(get(df), c("Vaccine","Dose","Delay","Day","2.5","25","50","75","97.5"))

VE <- rbind(VE_BB_prctile_7delay,VE_BB_prctile_14delay,VE_BB_prctile_21delay,VE_BBB_prctile_7delay,VE_BBB_prctile_14delay,VE_BBB_prctile_21delay,VE_BBBB_prctile_7delay,VE_BBBB_prctile_14delay,VE_BBBB_prctile_21delay,
            VE_CC_prctile_7delay,VE_CC_prctile_14delay,VE_CC_prctile_21delay,VE_CCC_prctile_7delay,VE_CCC_prctile_14delay,VE_CCC_prctile_21delay,VE_CCCC_prctile_7delay,VE_CCCC_prctile_14delay,VE_CCCC_prctile_21delay)
VE <- data.table::setnames(VE,c("Vaccine","Doses","Delay","Day","LCI","25","VE","75","UCI"))

VE$Delay <- fct_relevel(VE$Delay, "7 days","14 days","21 days")


VE$Vaccine <- as.factor(VE$Vaccine)
VE$Doses <- as.factor(VE$Doses)
VE$Delay <- as.factor(VE$Delay)


VE_plot <- VE %>% ggplot(aes(x=Day))+
  geom_line(aes(y=VE,colour=Delay),size=1)+
  geom_ribbon(aes(ymin=LCI,ymax=UCI,fill=Delay),alpha = 0.3)+  
  facet_grid(Doses~Vaccine)+
  theme_classic()+
  labs(x="Days from receiving last dose",y="Vaccine effectiveness",fill="Delay to vaccine effectiveness taking effect",colour="Delay to vaccine effectiveness taking effect")+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=(c(0,1)))+
  theme(legend.position="bottom",
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        panel.border=element_rect(colour='black',fill=NA,size=1),
        strip.text.x = element_text(size=8),
        strip.text.y = element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.spacing=unit(2,"lines")
        )
  
ggsave(file = "Extended_Data_Fig_4.tiff",VE_plot, height=5,width=5,dpi=1000,compression="lzw")


