library (tidyverse)
library (forcats)
library(cowplot)
library(lubridate)

rm(list=ls())
# dev.off()

## Purpose: To create Figure 3, 4, 5 and 6, and Extended data Figures 5, 6 and 7,
## and calculate ascertainment rate based on cumulative estimated and confirmed cases


setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
startDate = ymd("2022-01-01")
endDate = ymd("2022-07-31")


# VE plots (Fig 3)

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
        strip.text.x = element_text(size =7),
        strip.text.y = element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.spacing=unit(2,"lines")
  )


dev.new (width = 5, height = 7.5, unit = "px", noRStudioGD = TRUE)
png(file = "VE.png", width = 500, height = 750)

VE_plot

dev.off()


## Collect and wrangle data for IAR, Population Immunity and Ascertainment Rate plots

paths <- dir("./iar_imm_age/",pattern="\\.csv$",full.names = FALSE)
names(paths) <- basename(paths)
estimates <- map_dfr(paste0(getwd(),"/iar_imm_age/",paths),read.csv,stringsAsFactors=TRUE, header = FALSE)


names(estimates) <- c("mean","lci","Q1", "median", "Q3", "uci","Cohort","Delay","Measure")
parameter_levels <- c("0-11","12-19*","20-29", "30-39", "40-49", "50-59", "60+*", "all")
standard_weeks <- rep(c(seq (dmy("01/01/2022"),dmy("30/7/2022"),by="weeks"),dmy("31/07/2022")),length(paths))

estimates$Cohort <- as.factor(estimates$Cohort)
estimates$Cohort <- fct_relevel(estimates$Cohort,"0-11","12-19","20-29", "30-39", "40-49", "50-59", "60+", "all")
estimates$Delay <- as.factor(estimates$Delay)
estimates$Delay <- fct_relevel(estimates$Delay,"7 day delay","14 day delay","21 day delay")
estimates$Measure <- as.factor(estimates$Measure)
estimates$Measure <- fct_relevel(estimates$Measure,"IAR","Pop immunity","Pop immunity (Scenario 1)","Pop immunity (Scenario 2)")
estimates$Week <- standard_weeks

totalPopu <- read_csv("./VE/total_Popu.csv")
totalPopu <- totalPopu %>%
  group_by(cohort) %>%
  summarise(popu=sum(totalPopu)*1e3)
totalPopu <- rename(totalPopu,`Age group` = cohort)
totalPopu[nrow(totalPopu)+1,] <- list("all",sum(totalPopu[1:nrow(totalPopu),2]))
names(totalPopu) <- c("Cohort","Population")

popu_by_cohort <- inner_join(estimates,totalPopu,by="Cohort") %>% select(Population)

estimates_popu_level <- estimates
estimates_popu_level[,c(1:6)] <- sweep(as.matrix(estimates %>% select(mean,lci,Q1,median,Q3,uci)),MARGIN=1,as.vector(popu_by_cohort)[[1]],`*`)

IAR_popu_level <- estimates_popu_level %>% filter(Measure=="IAR")
IAR_popu_level_incre <- IAR_popu_level
IAR_popu_level_incre[2:nrow(IAR_popu_level_incre),c(1:6)] <-as.data.frame(pmax(as.matrix(IAR_popu_level[2:nrow(IAR_popu_level),c(1:6)]-IAR_popu_level[1:nrow(IAR_popu_level)-1,c(1:6)]),0))

total_cases <- read_csv("total_cases.csv")

estimated_cases_total <- estimates_popu_level %>% 
  filter(Measure == "IAR",Week == ymd("2022-07-31")) %>% 
  select(mean,lci,Q1,median,Q3,uci,Cohort,Delay)


ascertained_cases <- inner_join(estimated_cases_total,total_cases,by="Cohort")
ascertained_cases[,c(1:6)] <- sweep(1/as.matrix(ascertained_cases %>% select(mean,lci,Q1,median,Q3,uci)),MARGIN=1,
                                    as.vector(ascertained_cases %>% select(Cases))[[1]],`*`)
ascertained_cases$Cohort <- as.factor(ascertained_cases$Cohort)
ascertained_cases$test_type <- as.factor(ascertained_cases$test_type)

# IAR plot by age group (Fig 4 Panel A)
iar_cumu_plot <- estimates %>% filter(Measure == "IAR" & Week==ymd("2022-07-31")) %>%
  ggplot(aes(Cohort,median,colour=Delay)) +
  geom_pointrange (aes(y=median, ymin = lci, ymax = uci, colour = Delay),
                   position = position_dodge(0.5), size = 0.1)+
  labs(colour = "Assumption on delay to VE taking effect:",x= "Age group", y = "Infection attack rate")+
  scale_colour_manual(values = c("7 day delay" = "#E64B35FF", "14 day delay" = "#4DBBD5FF", "21 day delay" = "91D1C2FF"))+
  theme_classic()+
  theme(axis.title.y = element_text(size=15), axis.title.x = element_text(size=15), axis.text = element_text(size=15),
        legend.position = c(0.8,0.8), legend.title=element_text(size = 20), legend.text = element_text(size = 20)) +
  coord_cartesian(ylim=c(0.0,0.8))+
  scale_x_discrete (labels = parameter_levels)

iar_cumu_plot

# Population immunity plot by age group (Fig 4 Panel B) (no waning)
imm_cumu_plot <- estimates %>% filter(Measure == "Pop immunity" & Week==ymd("2022-07-31")) %>%
  ggplot(aes(Cohort,median,colour=Delay)) +
  geom_pointrange (aes(y=median, ymin = lci, ymax = uci, colour = Delay),
                   position = position_dodge(0.5), size = 0.1)+
  labs(colour = "Assumption on delay to VE taking effect:",x= "Age group", y = "Population immunity\n(no waning)")+
  scale_colour_manual(values = c("7 day delay" = "#E64B35FF", "14 day delay" = "#4DBBD5FF", "21 day delay" = "91D1C2FF"))+
  theme_classic()+
  theme(axis.title.y = element_text(size=15), axis.title.x = element_text(size=15), axis.text = element_text(size=15),
        legend.position = c(0.8,0.8), legend.title=element_text(size = 20), legend.text = element_text(size = 20)) +
  coord_cartesian(ylim=c(0.0,0.8))+
  scale_x_discrete (labels = parameter_levels)

imm_cumu_plot

# Ascertainment rate (overall) plot (Fig 4 Panel C)
ar_plot <- ascertained_cases %>% filter(test_type=="PCR+RAT") %>%
  ggplot(aes(Cohort,median,colour=Delay)) +
  geom_pointrange (aes(y=median, ymin = lci, ymax = uci, colour = Delay),
                   position = position_dodge(0.5), size = 0.1)+
  labs(colour = "Assumption on delay to VE taking effect:",x= "Age group", y = "Ascertainment rate")+
  scale_colour_manual(values = c("7 day delay" = "#E64B35FF", "14 day delay" = "#4DBBD5FF", "21 day delay" = "91D1C2FF"))+
  theme_classic()+
  theme(axis.title.y = element_text(size=15), axis.title.x = element_text(size=15), axis.text = element_text(size=15),
        legend.position = c(0.8,0.8), legend.title=element_text(size = 20), legend.text = element_text(size = 20)) +
  coord_cartesian(ylim=c(0.0,0.6))+
  scale_x_discrete (labels = parameter_levels)

ar_plot

legend_ar <- get_legend(
  ar_plot+
    guides(color = guide_legend(nrow=2, ncol =2),
           shape = guide_legend(override.aes = list(linetype = "blank")))+
    theme(legend.position = "bottom", legend.title=element_text(size=8), legend.text = element_text(size=8))
  # scale_colour_discrete(labels = c("No protection (VE = 0)","Same as BNT (VE = BNT VE)"))
)


ar_plot_with_legend <- plot_grid(ar_plot+theme(legend.position = "none"), legend_ar, ncol = 1, rel_heights = c(1,0.2))

ar_plot_with_legend

# Ascertainment rate (PCR) plot (Fig 2 Panel E)
ar_PCR_plot <- ascertained_cases %>% filter(test_type=="PCR") %>%
  ggplot(aes(Cohort,median,colour=Delay)) +
  geom_pointrange (aes(y=median, ymin = lci, ymax = uci, colour = Delay),
                   position = position_dodge(0.5), size = 0.1)+
  labs(colour = "Assumption on delay to VE taking effect:",x= "Age group", y = "Ascertainment rate \n(PCR only)")+
  scale_colour_manual(values = c("7 day delay" = "#E64B35FF", "14 day delay" = "#4DBBD5FF", "21 day delay" = "91D1C2FF"))+
  theme_classic()+
  theme(axis.title.y = element_text(size=15), axis.title.x = element_text(size=15), axis.text = element_text(size=15),
        legend.position = c(0.8,0.8), legend.title=element_text(size = 20), legend.text = element_text(size = 20)) +
  coord_cartesian(ylim=c(0.0,0.6))+
  scale_x_discrete (labels = parameter_levels)

ar_PCR_plot

legend_ar_pcr <- get_legend(
  ar_PCR_plot+
    guides(color = guide_legend(nrow=2, ncol =2),
           shape = guide_legend(override.aes = list(linetype = "blank")))+
    theme(legend.position = "bottom", legend.title=element_text(size=8), legend.text = element_text(size=8))
  # scale_colour_discrete(labels = c("No protection (VE = 0)","Same as BNT (VE = BNT VE)"))
)

ar_pcr_plot_with_legend <- plot_grid(ar_PCR_plot+theme(legend.position = "none"), legend_ar, ncol = 1, rel_heights = c(1,0.2))

ar_pcr_plot_with_legend

# IAR plot over time (all ages) (Figure 5 Panel A)
iar_plot <- estimates %>% filter(Cohort =="all"  & Measure == "IAR") %>% ggplot(aes(x=Week))+
  geom_line(aes(y=median,colour=Delay),size=1)+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Delay),alpha = 0.3)+  
  theme_classic()+
  labs(x="Date",y="Infection attack rate",fill="Delay to vaccine effectiveness \ntaking effect",colour="Delay to vaccine effectiveness \ntaking effect")+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=(c(0,0.75)))+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
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
iar_plot

# IAR plot over time (by age group) (Figure 6)
iar_plot_facet <- estimates %>% filter(Measure == "IAR") %>% ggplot(aes(x=Week))+
  geom_line(aes(y=median,colour=Delay),size=1)+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Delay),alpha = 0.3)+  
  theme_classic()+
  labs(x="Date",y="Infection attack rate",fill="Delay to vaccine effectiveness \ntaking effect",colour="Delay to vaccine effectiveness \ntaking effect")+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=(c(0,0.6)))+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
  theme(legend.position="bottom",
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        panel.border=element_rect(colour='black',fill=NA,size=1),
        strip.text.x = element_text(size=8),
        strip.text.y = element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.spacing=unit(2,"lines")
  )+
  facet_wrap(vars(Cohort),ncol=2)+
  theme(plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"))
iar_plot_facet

# Population immunity plot over time (all ages) (no waning) (Figure 5 Panel B)
imm_plot <- estimates %>% filter(Cohort == "all" & Measure == "Pop immunity") %>% ggplot(aes(x=Week))+
  geom_line(aes(y=median,colour=Delay),size=1)+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Delay),alpha = 0.3)+  
  theme_classic()+
  labs(x="Date",y="Population immunity\n(no waning)",fill="Delay to vaccine effectiveness \ntaking effect",colour="Delay to vaccine effectiveness \ntaking effect")+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=(c(0,0.75)))+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
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
imm_plot

# Population immunity plot over time (by age group) (no waning) (Figure 5 Panel C)
imm_plot_facet <- estimates %>% filter(Measure == "Pop immunity") %>% ggplot(aes(x=Week))+
  geom_line(aes(y=median,colour=Delay),size=1)+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Delay),alpha = 0.3)+  
  theme_classic()+
  labs(x="Date",y="Population immunity\n(no waning)",fill="Delay to vaccine effectiveness \ntaking effect",colour="Delay to vaccine effectiveness \ntaking effect")+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=(c(0,0.75)))+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
  theme(legend.position="bottom",
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        panel.border=element_rect(colour='black',fill=NA,size=1),
        strip.text.x = element_text(size=8),
        strip.text.y = element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.spacing=unit(2,"lines")
  )+
  facet_wrap(vars(Cohort),ncol=2)+
  theme(plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"))
imm_plot_facet

# Population immunity plot over time (all ages) (15% waning in 365 days) (Figure 5 Panel D) following [15] Barnard et al
imm_barnard_plot <- estimates %>% filter(Cohort == "all" & Measure == "Pop immunity (Scenario 1)") %>% ggplot(aes(x=Week))+
  geom_line(aes(y=median,colour=Delay),size=1)+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Delay),alpha = 0.3)+  
  theme_classic()+
  labs(x="Date",y="Population immunity\n(15% waning at 365 days)",fill="Delay to vaccine effectiveness \ntaking effect",colour="Delay to vaccine effectiveness \ntaking effect")+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=(c(0,0.75)))+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
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
imm_barnard_plot

# Population immunity plot over time (by age group) (15% waning in 365 days) (Extended Data Figure 5) following [15] Barnard et al
imm_barnard_plot_facet <- estimates %>% filter(Measure == "Pop immunity (Scenario 1)") %>% ggplot(aes(x=Week))+
  geom_line(aes(y=median,colour=Delay),size=1)+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Delay),alpha = 0.3)+  
  theme_classic()+
  labs(x="Date",y="Population immunity\n(15% waning at 365 days)",fill="Delay to vaccine effectiveness \ntaking effect",colour="Delay to vaccine effectiveness \ntaking effect")+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=(c(0,0.75)))+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
  theme(legend.position="bottom",
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        panel.border=element_rect(colour='black',fill=NA,size=1),
        strip.text.x = element_text(size=8),
        strip.text.y = element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.spacing=unit(2,"lines")
  )+
  facet_wrap(vars(Cohort),ncol=2)+
  theme(plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"))
imm_barnard_plot_facet

# Population immunity plot over time (all ages) (25% waning in 100 days) (Figure 5 Panel D) following [16] Malato et al and [17] Altarawneh et al
imm_malato_plot <- estimates %>% filter(Cohort == "all" & Measure == "Pop immunity (Scenario 2)") %>% ggplot(aes(x=Week))+
  geom_line(aes(y=median,colour=Delay),size=1)+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Delay),alpha = 0.3)+  
  theme_classic()+
  labs(x="Date",y="Population immunity\n(25% waning at 100 days)",fill="Delay to vaccine effectiveness \ntaking effect",colour="Delay to vaccine effectiveness \ntaking effect")+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=(c(0,0.75)))+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
  theme(legend.position="bottom",
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        panel.border=element_rect(colour='black',fill=NA,size=1),
        strip.text.x = element_text(size=8),
        strip.text.y = element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.spacing=unit(2,"lines")
  )+
  theme(plot.margin = margin(t = 5, r = 15, b = 5, l = 5, unit = "pt"))
imm_malato_plot

# Population immunity plot over time (by age group) (25% waning in 100 days) (Extended Data Figure 6) following [16] Malato et al and [17] Altarawneh et al
imm_malato_plot_facet <- estimates %>% filter(Measure == "Pop immunity (Scenario 2)") %>% ggplot(aes(x=Week))+
  geom_line(aes(y=median,colour=Delay),size=1)+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Delay),alpha = 0.3)+  
  theme_classic()+
  labs(x="Date",y="Population immunity\n(25% waning at 100 days)",fill="Delay to vaccine effectiveness \ntaking effect",colour="Delay to vaccine effectiveness \ntaking effect")+
  scale_fill_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_color_manual(values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=(c(0,0.75)))+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
  theme(legend.position="bottom",
        axis.title.x = element_text(size=8),
        axis.title.y = element_text(size=8),
        panel.border=element_rect(colour='black',fill=NA,size=1),
        strip.text.x = element_text(size=8),
        strip.text.y = element_text(size=8),
        legend.title = element_text(size=8),
        legend.text = element_text(size=8),
        panel.spacing=unit(2,"lines")
  )+
  facet_wrap(vars(Cohort),ncol=2)+
  theme(plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"))
imm_malato_plot_facet

iar_imm_plot <- align_plots(iar_plot+theme(legend.position="none",axis.text.x = element_text(angle=45,size=8,hjust=1),axis.text.y=element_text(size=8),axis.title.x = element_blank(),axis.title.y = element_text(size=8)),
                            imm_plot+theme(legend.position="none",axis.text.x = element_text(angle=45,size=8,hjust=1),axis.text.y=element_text(size=8),axis.title.x = element_blank(),axis.title.y = element_text(size=8)),
                            imm_barnard_plot+theme(legend.position="none",axis.text.x = element_text(angle=45,size=8,hjust=1),axis.text.y=element_text(size=8),axis.title.x = element_blank(),axis.title.y = element_text(size=8)),
                            imm_malato_plot+theme(legend.position="none",axis.text.x = element_text(angle=45,size=8,hjust=1),axis.text.y=element_text(size=8),axis.title.x = element_blank(),axis.title.y = element_text(size=8)),
                            align="h",axis="lr",greedy=TRUE)

iar_imm_plot_with_labels <- plot_grid(iar_imm_plot[[1]],iar_imm_plot[[2]],iar_imm_plot[[3]],iar_imm_plot[[4]],nrow=1,labels=c("A","B","C","D"))

iar_imm_cumu_plot <- align_plots(iar_cumu_plot+theme(legend.position="none",axis.text.x = element_text(size=8),axis.text.y=element_text(size=8),axis.title.x = element_blank(),axis.title.y = element_text(size=8)),
                                 imm_cumu_plot+theme(legend.position="none",axis.text.x = element_text(size=8),axis.text.y=element_text(size=8),axis.title.x = element_text(size=8),axis.title.y = element_text(size=8)),
                                 align="v",axis="lr")

iar_imm_cumu_plot_with_labels <- plot_grid(iar_imm_cumu_plot[[1]],iar_imm_cumu_plot[[2]],ncol=1,labels=c("A","B"),label_y = c(1,1))

ar_all_plots <- align_plots(ar_plot+theme(legend.position="none", axis.text.x = element_text(size=8),axis.text.y=element_text(size=8),axis.title.x = element_blank(),axis.title.y = element_text(size=8)),
                            ar_PCR_plot+theme(legend.position="none", axis.text.x = element_text(size=8),axis.text.y=element_text(size=8),axis.title.x = element_text(size=8),axis.title.y = element_text(size=8)),
                            align="v",axis="lr")
ar_all_plots_with_labels <- plot_grid(ar_all_plots[[1]],ar_all_plots[[2]],ncol=1,labels=c("C","D"),label_y=c(1, 1))


legend_VE <- get_legend(
  VE_plot+
    guides(color = guide_legend(nrow=1, ncol =1),
           shape = guide_legend(override.aes = list(linetype = "blank")))+
    theme(legend.position = "bottom", legend.title=element_text(size=8), legend.text = element_text(size = 8))
)

iar_imm_cumu_ar_plot <- plot_grid(iar_imm_cumu_plot_with_labels,ar_all_plots_with_labels)

full_plot_fig_3 <- plot_grid(
  VE_plot+theme(axis.text.x = element_text(size=8),axis.text.y=element_text(size=8),axis.title.x = element_text(size=8),axis.title.y = element_text(size=8))
)

full_plot_fig_4 <- plot_grid(
  iar_imm_cumu_ar_plot+theme(axis.title.x = element_text(size=8),axis.title.y = element_text(size=8)),
  legend_VE,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(1, 0.1))

full_plot_fig_5 <- plot_grid(
  iar_imm_plot_with_labels+theme(axis.text.x=element_blank(),axis.title.x = element_text(size=8),axis.title.y = element_text(size=8)),
  legend_VE,
  ncol = 1,
  # rel_heights values control vertical title margins
  rel_heights = c(1, 0.1))

ggsave(file = "Fig_3.pdf",full_plot_fig_3,width=120,height=120,units="mm")
ggsave(file = "Fig_4.pdf",full_plot_fig_4,width=180,height=120,units="mm")
ggsave(file = "Fig_5.pdf",full_plot_fig_5,width=180,height=80,units="mm")
ggsave(file = "Fig_6.pdf",iar_plot_facet,width=180,height=200,units="mm")

ggsave(file = "Extended_Data_Fig_5.tiff",imm_plot_facet,height=9,width=7,dpi=1000,compression="lzw")
ggsave(file = "Extended_Data_Fig_6.tiff",imm_barnard_plot_facet,height=9,width=7,dpi=1000,compression="lzw")
ggsave(file = "Extended_Data_Fig_7.tiff",imm_malato_plot_facet,height=9,width=7,dpi=1000,compression="lzw")

