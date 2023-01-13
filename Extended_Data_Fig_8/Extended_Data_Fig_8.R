library (tidyverse)
library (forcats)
library(cowplot)
library(lubridate)

rm(list=ls())

## Purpose: To create Extended Data Figure 8

setwd(dirname(rstudioapi::getActiveDocumentContext()$path))
startDate = ymd("2022-01-01")
endDate = ymd("2022-07-31")

# dev.off()

## Collect and wrangle data for IAR plots under each seroconversion rate assumption (100% (Base Case), 90% and 75%)

paths <- dir("./iar_age/",pattern="\\.csv$",full.names = FALSE)
names(paths) <- basename(paths)
estimates <- map_dfr(paste0(getwd(),"/iar_age/",paths),read.csv,stringsAsFactors=TRUE, header = FALSE)


names(estimates) <- c("mean","lci","Q1", "median", "Q3", "uci","Cohort","Delay","Measure")
parameter_levels <- c("0-11","12-19*","20-29", "30-39", "40-49", "50-59", "60+*", "all")
standard_weeks <- rep(c(seq (dmy("01/01/2022"),dmy("30/7/2022"),by="weeks"),dmy("31/07/2022")),length(paths))

estimates$Cohort <- as.factor(estimates$Cohort)
estimates$Cohort <- fct_relevel(estimates$Cohort,"0-11","12-19","20-29", "30-39", "40-49", "50-59", "60+", "all")
estimates$Delay <- as.factor(estimates$Delay)
estimates$Delay <- fct_relevel(estimates$Delay,"7 day delay")
estimates$Measure <- as.factor(estimates$Measure)
estimates$Measure <- fct_relevel(estimates$Measure,"IAR","IAR10","IAR25")
estimates$Week <- standard_weeks

# IAR plot over time (by age group and assumption on seroconversion rate) (Extended Data Figure 10)

iar_plot_facet <- estimates %>% ggplot(aes(x=Week))+
  geom_line(aes(y=median,colour=Measure),size=1)+
  geom_ribbon(aes(ymin=lci,ymax=uci,fill=Measure),alpha = 0.3)+  
  theme_classic()+
  labs(x="Date",y="Infection attack rate\n(7 day delay to VE taking effect)",fill="Assumption on seroconversion rate",colour="Assumption on seroconversion rate")+
  scale_fill_manual(labels=c("100% (Base Case)","90%","75%"),values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_color_manual(labels=c("100% (Base Case)","90%","75%"),values=c("#E64B35FF","#4DBBD5FF","#91D1C2FF"))+
  scale_y_continuous(expand=c(0,0))+
  coord_cartesian(ylim=(c(0,0.75)))+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
  theme(legend.position="bottom",
        axis.title.x = element_text(size=12),
        axis.title.y = element_text(size=12),
        panel.border=element_rect(colour='black',fill=NA,size=1),
        strip.text.x = element_text(size = 12),
        strip.text.y = element_text(size=12),
        legend.title = element_text(size=10),
        legend.text = element_text(size=10),
        panel.spacing=unit(2,"lines")
  )+
  facet_wrap(vars(Cohort),ncol=2)+
  theme(plot.margin = margin(t = 5, r = 20, b = 5, l = 5, unit = "pt"))
iar_plot_facet

ggsave(file = "Extended_Data_Fig_8.tiff",iar_plot_facet,height=9,width=7,dpi=1000,compression="lzw")

