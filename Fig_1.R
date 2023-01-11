library(tidyverse)
library(lubridate)
library(cowplot)
library(janitor)

rm(list=ls())
setwd(dirname(rstudioapi::getActiveDocumentContext()$path))

# Purpose: Plot Figure 1
# Note 1: Vaccination history and viral load data from sewage surveillance cannot be shared due to 
# confidentiality undertaking with the Department of Health and the Environmental Protection Department, both
# of the Government of the HKSAR. 

# Plot vaccination histories (Panels D-K)
vax_startDate = ymd("2021-01-01")
startDate = ymd("2022-01-01")
endDate = ymd("2022-07-31")

BNT_date <- read_csv("combined_data.csv")
BNT_date <- BNT_date %>% 
  filter(no_of_vaccines>=1 & vaccine_category == "BNT") %>%
  select(no_of_vaccines, dose_1_vaccination_date,dose_2_vaccination_date,dose_3_vaccination_date,dose_4_vaccination_date)
BNT_date <- BNT_date[rowSums(is.na(BNT_date))!=ncol(BNT_date),]
colnames(BNT_date) <- c("no_of_vaccines", "dose_1","dose_2","dose_3","dose_4")
BNT_date$dose_1 <- as.Date(BNT_date$dose_1)
BNT_date$dose_2 <- as.Date(BNT_date$dose_2)
BNT_date$dose_3 <- as.Date(BNT_date$dose_3)
BNT_date$dose_4 <- as.Date(BNT_date$dose_4)

BNT_date_dose_1 <- BNT_date %>% filter(no_of_vaccines == 1) %>% select(dose_1)
BNT_date_dose_2 <- BNT_date %>% filter(no_of_vaccines == 2) %>% select(dose_1,dose_2)
BNT_date_dose_3 <- BNT_date %>% filter(no_of_vaccines == 3) %>% select(dose_1,dose_2,dose_3)
BNT_date_dose_4 <- BNT_date %>% filter(no_of_vaccines >= 4) %>% select(dose_1,dose_2,dose_3,dose_4)

# Dose 1 BNT 
BNT_date_dose_1 <- BNT_date_dose_1 %>% pivot_longer(everything(),names_to = "Dose",values_to = "Date")
BNT_count_dose_1 <- count(BNT_date_dose_1,Date,Dose)
BNT_count_dose_1 <- BNT_count_dose_1[!is.na(BNT_count_dose_1$Date) & BNT_count_dose_1$Date >= vax_startDate & BNT_count_dose_1$Date <= endDate,]
BNT_count_by_week_dose_1 <- BNT_count_dose_1 %>% mutate(year_week = floor_date(Date,"1 week")) %>%
  group_by (year_week,Dose) %>%
  summarize (weekly_total = sum(n))

BNT_count_by_week_dose_1 <- BNT_count_by_week_dose_1 %>%
  mutate(Dose=
           case_when(
             Dose=="dose_1" ~ "Dose 1",
           ))

BNT_barplot_dose_1 <- BNT_count_by_week_dose_1 %>% ggplot(aes(x=year_week, y = weekly_total, fill=Dose)) + 
  annotate(geom = "rect",xmin=startDate,xmax=endDate,ymin=0,ymax=Inf,fill="#F39B7F99",alpha=0.1)+
  geom_col() + 
  labs(y = "Subjects")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months", limits = as.Date(c(vax_startDate,endDate)))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2"))+
  scale_y_continuous(limits=c(0,15),expand=c(0,0))+
  scale_fill_manual(values=c('Dose 1' = "#E64B35FF"))+
  theme(axis.title.x = element_blank())+
  theme_classic()+
  annotate("text",x=dmy("15/04/2022"),y=12,label="Fifth wave",size=2.5)

BNT_barplot_dose_1

# Dose 2 BNT 
BNT_date_dose_2 <- BNT_date_dose_2 %>% pivot_longer(everything(),names_to = "Dose",values_to = "Date")
BNT_count_dose_2 <- count(BNT_date_dose_2,Date,Dose)
BNT_count_dose_2 <- BNT_count_dose_2[!is.na(BNT_count_dose_2$Date) & BNT_count_dose_2$Date >= vax_startDate & BNT_count_dose_2$Date <= endDate,]
BNT_count_by_week_dose_2 <- BNT_count_dose_2 %>% mutate(year_week = floor_date(Date,"1 week")) %>%
  group_by (year_week,Dose) %>%
  summarize (weekly_total = sum(n))

BNT_count_by_week_dose_2 <- BNT_count_by_week_dose_2 %>%
  mutate(Dose=
           case_when(
             Dose=="dose_1" ~ "Dose 1",
             Dose=="dose_2" ~ "Dose 2"
           ))

BNT_barplot_dose_2 <- BNT_count_by_week_dose_2 %>% ggplot(aes(x=year_week, y = weekly_total, fill=Dose)) + 
  annotate(geom = "rect",xmin=startDate,xmax=endDate,ymin=0,ymax=Inf,fill="#F39B7F99",alpha=0.1)+
  geom_col() + 
  labs(y = "Subjects")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months", limits = as.Date(c(vax_startDate,endDate)))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2"))+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))+
  scale_fill_manual(values=c('Dose 1' = "#E64B35FF", 'Dose 2' = "#4DBBD5FF"))+
  theme(axis.title.x = element_blank())+
  theme_classic()

BNT_barplot_dose_2

# Dose 3 BNT
BNT_date_dose_3 <- BNT_date_dose_3 %>% pivot_longer(everything(),names_to = "Dose",values_to = "Date")
BNT_count_dose_3 <- count(BNT_date_dose_3,Date,Dose)
BNT_count_dose_3 <- BNT_count_dose_3[!is.na(BNT_count_dose_3$Date) & BNT_count_dose_3$Date >= vax_startDate & BNT_count_dose_3$Date <= endDate,]
BNT_count_by_week_dose_3 <- BNT_count_dose_3 %>% mutate(year_week = floor_date(Date,"1 week")) %>%
  group_by (year_week,Dose) %>%
  summarize (weekly_total = sum(n))

BNT_count_by_week_dose_3 <- BNT_count_by_week_dose_3 %>%
  mutate(Dose=
           case_when(
             Dose=="dose_1" ~ "Dose 1",
             Dose=="dose_2" ~ "Dose 2",
             Dose=="dose_3" ~ "Dose 3"
           ))

BNT_barplot_dose_3 <- BNT_count_by_week_dose_3 %>% ggplot(aes(x=year_week, y = weekly_total, fill=Dose)) + 
  annotate(geom = "rect",xmin=startDate,xmax=endDate,ymin=0,ymax=Inf,fill="#F39B7F99",alpha=0.1)+
  geom_col() + 
  labs(y = "Subjects")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months", limits = as.Date(c(vax_startDate,endDate)))+
  scale_y_continuous(limits=c(0,300),expand=c(0,0))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2","Dose 3"))+
  scale_fill_manual(values = c("Dose 1" = "#E64B35FF", "Dose 2" = "#4DBBD5FF",
                               "Dose 3" = "#00A087FF"))+
  theme_classic()

BNT_barplot_dose_3

# Dose 4 BNT
BNT_date_dose_4 <- BNT_date_dose_4 %>% pivot_longer(everything(),names_to = "Dose",values_to = "Date")
BNT_count_dose_4 <- count(BNT_date_dose_4,Date,Dose)
BNT_count_dose_4 <- BNT_count_dose_4[!is.na(BNT_count_dose_4$Date) & BNT_count_dose_4$Date >= vax_startDate & BNT_count_dose_4$Date <= endDate,]
BNT_count_by_week_dose_4 <- BNT_count_dose_4 %>% mutate(year_week = floor_date(Date,"1 week")) %>%
  group_by (year_week,Dose) %>%
  summarize (weekly_total = sum(n))

BNT_count_by_week_dose_4 <- BNT_count_by_week_dose_4 %>%
  mutate(Dose=
           case_when(
             Dose=="dose_1" ~ "Dose 1",
             Dose=="dose_2" ~ "Dose 2",
             Dose=="dose_3" ~ "Dose 3",
             Dose=="dose_4" ~ "Dose 4"
           ))

BNT_barplot_dose_4 <- BNT_count_by_week_dose_4 %>% ggplot(aes(x=year_week, y = weekly_total, fill=Dose)) + 
  annotate(geom = "rect",xmin=startDate,xmax=endDate,ymin=0,ymax=Inf,fill="#F39B7F99",alpha=0.1)+
  geom_col() + 
  labs(x = "Date", y = "Subjects")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months", limits = as.Date(c(vax_startDate,endDate)))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2","Dose 3", "Dose 4"))+
  scale_y_continuous(limits=c(0,15),expand=c(0,0))+
  scale_fill_manual(values = c("Dose 1" = "#E64B35FF", "Dose 2" = "#4DBBD5FF",
                               "Dose 3" = "#00A087FF", "Dose 4" = "#3C5488FF"))+
  theme_classic()
BNT_barplot_dose_4

C_date <- read_csv("combined_data.csv")
C_date <- C_date %>% 
  filter(no_of_vaccines>=1 & vaccine_category == "Sinovac") %>%
  select(no_of_vaccines, dose_1_vaccination_date,dose_2_vaccination_date,dose_3_vaccination_date,dose_4_vaccination_date)
C_date <- C_date[rowSums(is.na(C_date))!=ncol(C_date),]
colnames(C_date) <- c("no_of_vaccines", "dose_1","dose_2","dose_3","dose_4")
C_date$dose_1 <- as.Date(C_date$dose_1)
C_date$dose_2 <- as.Date(C_date$dose_2)
C_date$dose_3 <- as.Date(C_date$dose_3)
C_date$dose_4 <- as.Date(C_date$dose_4)

C_date_dose_1 <- C_date %>% filter(no_of_vaccines == 1) %>% select(dose_1)
C_date_dose_2 <- C_date %>% filter(no_of_vaccines == 2) %>% select(dose_1,dose_2)
C_date_dose_3 <- C_date %>% filter(no_of_vaccines == 3) %>% select(dose_1,dose_2,dose_3)
C_date_dose_4 <- C_date %>% filter(no_of_vaccines >= 4) %>% select(dose_1,dose_2,dose_3,dose_4)

# Dose 1 CoronaVac
C_date_dose_1 <- C_date_dose_1 %>% pivot_longer(everything(),names_to = "Dose",values_to = "Date")
C_count_dose_1 <- count(C_date_dose_1,Date,Dose)
C_count_dose_1 <- C_count_dose_1[!is.na(C_count_dose_1$Date) & C_count_dose_1$Date >= vax_startDate & C_count_dose_1$Date <= endDate,]
C_count_by_week_dose_1 <- C_count_dose_1 %>% mutate(year_week = floor_date(Date,"1 week")) %>%
  group_by (year_week,Dose) %>%
  summarize (weekly_total = sum(n))

C_count_by_week_dose_1 <- C_count_by_week_dose_1 %>%
  mutate(Dose=
           case_when(
             Dose=="dose_1" ~ "Dose 1",
           ))

C_barplot_dose_1 <- C_count_by_week_dose_1 %>% ggplot(aes(x=year_week, y = weekly_total, fill=Dose)) + 
  annotate(geom = "rect",xmin=startDate,xmax=endDate,ymin=0,ymax=Inf,fill="#F39B7F99",alpha=0.1)+
  geom_col() + 
  labs(y = "Subjects")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months", limits = as.Date(c(vax_startDate,endDate)))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2"))+
  scale_y_continuous(limits=c(0,15),expand=c(0,0))+
  scale_fill_manual(values=c('Dose 1' = "#E64B35FF"))+
  theme(axis.title.x = element_blank())+
  theme_classic()+
  annotate("text",x=dmy("15/04/2022"),y=12,label="Fifth wave",size=2.5)

C_barplot_dose_1

# Dose 2 CoronaVac
C_date_dose_2 <- C_date_dose_2 %>% pivot_longer(everything(),names_to = "Dose",values_to = "Date")
C_count_dose_2 <- count(C_date_dose_2,Date,Dose)
C_count_dose_2 <- C_count_dose_2[!is.na(C_count_dose_2$Date) & C_count_dose_2$Date >= vax_startDate & C_count_dose_2$Date <= endDate,]
C_count_by_week_dose_2 <- C_count_dose_2 %>% mutate(year_week = floor_date(Date,"1 week")) %>%
  group_by (year_week,Dose) %>%
  summarize (weekly_total = sum(n))

C_count_by_week_dose_2 <- C_count_by_week_dose_2 %>%
  mutate(Dose=
           case_when(
             Dose=="dose_1" ~ "Dose 1",
             Dose=="dose_2" ~ "Dose 2"
           ))

C_barplot_dose_2 <- C_count_by_week_dose_2 %>% ggplot(aes(x=year_week, y = weekly_total, fill=Dose)) + 
  annotate(geom = "rect",xmin=startDate,xmax=endDate,ymin=0,ymax=Inf,fill="#F39B7F99",alpha=0.1)+
  geom_col() + 
  labs(y = "Subjects")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months", limits = as.Date(c(vax_startDate,endDate)))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2"))+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))+
  scale_fill_manual(values=c('Dose 1' = "#E64B35FF", 'Dose 2' = "#4DBBD5FF"))+
  theme(axis.title.x = element_blank())+
  theme_classic()

C_barplot_dose_2

# Dose 3 CoronaVac
C_date_dose_3 <- C_date_dose_3 %>% pivot_longer(everything(),names_to = "Dose",values_to = "Date")
C_count_dose_3 <- count(C_date_dose_3,Date,Dose)
C_count_dose_3 <- C_count_dose_3[!is.na(C_count_dose_3$Date) & C_count_dose_3$Date >= vax_startDate & C_count_dose_3$Date <= endDate,]
C_count_by_week_dose_3 <- C_count_dose_3 %>% mutate(year_week = floor_date(Date,"1 week")) %>%
  group_by (year_week,Dose) %>%
  summarize (weekly_total = sum(n))

C_count_by_week_dose_3 <- C_count_by_week_dose_3 %>%
  mutate(Dose=
           case_when(
             Dose=="dose_1" ~ "Dose 1",
             Dose=="dose_2" ~ "Dose 2",
             Dose=="dose_3" ~ "Dose 3"
           ))

C_barplot_dose_3 <- C_count_by_week_dose_3 %>% ggplot(aes(x=year_week, y = weekly_total, fill=Dose)) + 
  annotate(geom = "rect",xmin=startDate,xmax=endDate,ymin=0,ymax=Inf,fill="#F39B7F99",alpha=0.1)+
  geom_col() + 
  labs(y = "Subjects")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months", limits = as.Date(c(vax_startDate,endDate)))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2","Dose 3"))+
  scale_y_continuous(limits=c(0,100),expand=c(0,0))+
  scale_fill_manual(values = c("Dose 1" = "#E64B35FF", "Dose 2" = "#4DBBD5FF",
                               "Dose 3" = "#00A087FF"))+
  theme_classic()

C_barplot_dose_3

# Dose 4 CoronaVac
C_date_dose_4 <- C_date_dose_4 %>% pivot_longer(everything(),names_to = "Dose",values_to = "Date")
C_count_dose_4 <- count(C_date_dose_4,Date,Dose)
C_count_dose_4 <- C_count_dose_4[!is.na(C_count_dose_4$Date) & C_count_dose_4$Date >= vax_startDate & C_count_dose_4$Date <= endDate,]
C_count_by_week_dose_4 <- C_count_dose_4 %>% mutate(year_week = floor_date(Date,"1 week")) %>%
  group_by (year_week,Dose) %>%
  summarize (weekly_total = sum(n))

C_count_by_week_dose_4 <- C_count_by_week_dose_4 %>%
  mutate(Dose=
           case_when(
             Dose=="dose_1" ~ "Dose 1",
             Dose=="dose_2" ~ "Dose 2",
             Dose=="dose_3" ~ "Dose 3",
             Dose=="dose_4" ~ "Dose 4"
           ))

C_barplot_dose_4 <- C_count_by_week_dose_4 %>% ggplot(aes(x=year_week, y = weekly_total, fill=Dose)) + 
  annotate(geom = "rect",xmin=startDate,xmax=endDate,ymin=0,ymax=Inf,fill="#F39B7F99",alpha=0.1)+
  geom_col() + 
  labs(x = "Date", y = "Subjects")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "3 months", limits = as.Date(c(vax_startDate,endDate)))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2","Dose 3", "Dose 4"))+
  scale_y_continuous(limits=c(0,15),expand=c(0,0))+
  scale_fill_manual(values = c("Dose 1" = "#E64B35FF", "Dose 2" = "#4DBBD5FF",
                               "Dose 3" = "#00A087FF", "Dose 4" = "#3C5488FF"))+
  theme_classic()
C_barplot_dose_4

# Plot viral load vs incidence (Panel A)

viral_load <- read_csv("viral_load.csv")%>% 
  clean_names()

viral_load[viral_load$date==dmy("2/7/2022"),'x2_day_running_geometric_mean_viral_load'] <- viral_load[viral_load$date==dmy("1/7/2022"),'x2_day_running_geometric_mean_viral_load']

viral_load <- viral_load %>%
  select(date, x2_day_running_geometric_mean_viral_load) %>%
  filter(!is.na(x2_day_running_geometric_mean_viral_load))

colnames(viral_load) <- c("date","two_day_viral_load")


viral_load <- viral_load %>%
  mutate(proxy_day=as.numeric(date-ymd(startDate)+1,"days"))


start_viral_load <- data.frame(seq(startDate,ymd("2022-01-01"),by='days'),rep(0,ymd("2022-01-01")-startDate+1),1)
colnames(start_viral_load) <-c("date","two_day_viral_load","proxy_day")
viral_load <- rbind(start_viral_load,viral_load[,c('date','two_day_viral_load','proxy_day')])
viral_load$two_day_viral_load = viral_load$two_day_viral_load/1000000
viral_load <- viral_load %>% filter(date>=startDate & date<=endDate)

total_confirmed_cases <- read_csv("confirmed_cases.csv")
total_confirmed_cases$PCR <- rowSums(cbind(total_confirmed_cases$`Number of confirmed cases`,
                                                   total_confirmed_cases$`Number of cases tested positive for SARS-CoV-2 virus by nucleic acid tests`),na.rm = TRUE)
total_confirmed_cases$RAT <- rowSums(cbind(rep(0,nrow(total_confirmed_cases)),total_confirmed_cases$`Number of cases tested positive for SARS-CoV-2 virus by rapid antigen tests`),na.rm=TRUE)
                                     
confirmed_cases_dates <- total_confirmed_cases[2:nrow(total_confirmed_cases),1]
confirmed_daily_PCR <- total_confirmed_cases$PCR[2:nrow(total_confirmed_cases)]-total_confirmed_cases$PCR[1:nrow(total_confirmed_cases)-1]
confirmed_daily_RAT <- total_confirmed_cases$RAT[2:nrow(total_confirmed_cases)]-total_confirmed_cases$RAT[1:nrow(total_confirmed_cases)-1]

confirmed_cases <- cbind(confirmed_cases_dates,confirmed_daily_PCR,confirmed_daily_RAT)
names(confirmed_cases) <- c("Date","Daily confirmed cases from PCR","Daily confirmed cases from RAT")
confirmed_cases$Date <- dmy(confirmed_cases$Date)
confirmed_cases <- confirmed_cases %>% filter(Date >= startDate & Date <= endDate)
# confirmed_cases <- confirmed_cases %>% pivot_longer(!Date,names_to="Test type",values_to = "Daily Cases")
  
viral_load$`Daily confirmed\ncases (RAT)` <- unlist(confirmed_cases %>% select(`Daily confirmed cases from RAT`))
viral_load$`Daily confirmed\ncases (PCR)` <- unlist(confirmed_cases %>% select(`Daily confirmed cases from PCR`))
viral_load <- viral_load %>% pivot_longer(cols=starts_with("Daily"),names_to="Test type",values_to="daily_cases")
viral_load$`Test type` <- as.factor(viral_load$`Test type`)
viral_load <- viral_load %>% mutate (`Test type` = fct_relevel(`Test type`,"Daily confirmed\ncases (RAT)","Daily confirmed\ncases (PCR)"))

# viral_load$two_day_viral_load = viral_load$two_day_viral_load/sum(viral_load$two_day_viral_load)
coeff=1e4

viral_load_lineplot <- viral_load %>% ggplot(aes(x=date))+
  geom_col(aes(y=daily_cases/coeff,fill=`Test type`),width=1)+
  geom_line(aes(y=two_day_viral_load,colour="Viral load\nfrom wastewater"))+
  labs(y=expression(paste(10^6," RNA copies/L")),colour=NULL,fill=NULL)+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
  scale_y_continuous(expand=c(0,0),name=expression(paste(10^6, "RNA copies/L")),
                     sec.axis = sec_axis(~.*coeff,name="Confirmed cases per day"))+
  coord_cartesian(ylim = c(0, 8)) +
  scale_fill_manual(values=c("#91D1C2FF","#3C5488FF"))+
  scale_color_manual(values="#DC0000FF")+
  theme_classic()+
  theme(axis.title.x=element_blank(),legend.position=c(0.7,0.5),legend.text =element_text(size=7),legend.spacing.y = unit(-0.1,"cm"))+
  guides(fill=guide_legend(byrow = TRUE))

viral_load_lineplot

# Plot weekly sequences by lineage (from GISAID data) (Panel B)

sequences <- read_csv("SeqhkLineage.csv")
sequences <- sequences %>% select(date,lineage)
colnames(sequences) <- c("date","Lineage")
sequences <- sequences %>% mutate (Lineage = ifelse(Lineage=="B.1.1.529"|Lineage=="BA.3","Other",Lineage))
sequences_count <- count(sequences,date,Lineage)

sequences_count <- sequences_count %>% mutate(year_week=floor_date(date,"1 week")) %>%
  group_by(year_week,Lineage) %>%
  summarise(weekly_total=sum(n))

sequences_count <-sequences_count %>%
  group_by(year_week) %>%
  mutate(per=prop.table(weekly_total))

sequences_count_plot <- sequences_count %>% ggplot(aes(x=year_week, y = per, fill=Lineage)) + 
  geom_col() + 
  labs(x = "Date", y = "Proportion")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2","Dose 3", "Dose 4"))+
  scale_y_continuous(limits=c(0,1),expand=c(0,0))+
  scale_fill_manual(values = c("Delta" = "#E64B35FF", "BA.1" = "#4DBBD5FF",
                               "BA.2" = "#00A087FF", "BA.4" = "91D1C2FF",
                                 "BA.5" = "#F39B7FFF","Other" = "#3C5488FF"))+
  theme_classic()
sequences_count_plot

# Plot weekly cumulative vaccinations (by dose) (Panel C)
vaccines <- read_csv ("hk_vaccines.csv")
vaccines <- vaccines %>% select(date,firstDose.cumulative.total,secondDose.cumulative.total,thirdDose.cumulative.total,fourthDose.cumulative.total) 
vaccinees_by_dose <- as.data.frame(vaccines$date)
vaccinees_by_dose <- rename(vaccinees_by_dose,date = `vaccines$date`)

  
vaccinees_by_dose$`Dose 1` <- vaccines$firstDose.cumulative.total-vaccines$secondDose.cumulative.total
vaccinees_by_dose$`Dose 2` <- vaccines$secondDose.cumulative.total-vaccines$thirdDose.cumulative.total
vaccinees_by_dose$`Dose 3` <- vaccines$thirdDose.cumulative.total-vaccines$fourthDose.cumulative.total
vaccinees_by_dose$`Dose 4` <- vaccines$fourthDose.cumulative.total

vaccinees_by_dose <- vaccinees_by_dose %>% pivot_longer(!date,names_to = "Dose",values_to = "Count")
vaccinees_by_dose$Dose <- as.factor(vaccinees_by_dose$Dose)
vaccinees_by_dose <- vaccinees_by_dose %>% mutate (Dose = fct_relevel(Dose,"Dose 4","Dose 3","Dose 2","Dose 1"))
vaccinees_by_dose <- vaccinees_by_dose %>% filter (date >= startDate & date <= endDate )

vaccinees_by_dose_plot <- vaccinees_by_dose %>% ggplot(aes(x=date, y = Count/10^6, fill=Dose)) + 
  geom_col() + 
  labs(x = "Date", y = "Population\n(millions)",width=1)+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
  # scale_fill_discrete(labels = c("Dose 1", "Dose 2","Dose 3", "Dose 4"))+
  scale_y_continuous(limits=c(0,8),expand=c(0,0))+
  scale_fill_manual(values = c("Dose 1" = "#E64B35FF", "Dose 2" = "#4DBBD5FF",
                               "Dose 3" = "#00A087FF", "Dose 4" = "#3C5488FF"))+
  theme_classic()
  
vaccinees_by_dose_plot

# Plot weekly serum samples collected (Panel D)
samples <- read_csv("combined_data.csv") %>% select(donation_date)
samples_count <- count(samples,donation_date)

samples_count <- samples_count %>% mutate(year_week=floor_date(donation_date,"1 week")) %>%
  group_by(year_week) %>%
  summarise(weekly_total=sum(n))

samples_count_plot <- samples_count %>% ggplot(aes(x=year_week, y = weekly_total)) + 
  geom_col(fill="#3C548899") + 
  labs(x = "Date", y = "Weekly\nsamples")+
  scale_x_date(date_labels = "%Y-%m", date_breaks = "2 months", limits = as.Date(c(startDate,endDate)))+
  scale_y_continuous(limits=c(0,1000),expand=c(0,0))+
  theme_classic()
samples_count_plot

# Combine all plots and output as Fig_1.pdf

legend <- get_legend(
  BNT_barplot_dose_4+
    guides(colour=guide_legend(nrow=1))+
    theme(legend.position = "bottom",legend.box.margin = margin(0,0,0,0),legend.title = element_text(size=7), legend.text = element_text(size=7))
)


aligned_incidence_plot <- align_plots(viral_load_lineplot+theme(axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7),legend.title = element_text(size=7), legend.text = element_text(size=7)),
                                      sequences_count_plot+theme(axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7),legend.title = element_text(size=7), legend.text = element_text(size=7)),
                                      vaccinees_by_dose_plot+theme(axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7),legend.title = element_text(size=7), legend.text = element_text(size=7)),
                                      samples_count_plot+theme(axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7),legend.title = element_text(size=7), legend.text = element_text(size=7)),
                                      align="v",axis="lr")
dev.new(width=5,height=10,unit="px",noRStudioGD = TRUE)

png(file = "Figure_1a.png", width = 600, height = 400)
incidence_variant_plot <- plot_grid(aligned_incidence_plot[[1]],aligned_incidence_plot[[2]],aligned_incidence_plot[[3]],aligned_incidence_plot[[4]],ncol=2,nrow=2,rel_widths =c(1,1),rel_heights=c(1,1),labels=c("A","B","C","D"),label_y=c(1,1,1,1),label_size=10)
incidence_variant_plot
dev.off()

aligned_plots_BNT <- align_plots((BNT_barplot_dose_1+theme(legend.position="none",axis.text.x = element_text(size=7),axis.text.y=element_text(size=7), axis.title.x = element_blank(),axis.title.y = element_text(size=7))),
                             (BNT_barplot_dose_2+theme(legend.position="none",axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7))),
                             (BNT_barplot_dose_3+theme(legend.position="none",axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7))),
                             (BNT_barplot_dose_4+theme(legend.position="none",axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7))),
                             align="v",axis="lr")

vax_plots_BNT <- plot_grid(aligned_plots_BNT[[1]],aligned_plots_BNT[[2]],aligned_plots_BNT[[3]],aligned_plots_BNT[[4]],
          ncol=1, rel_heights = c(1,1,1,1),labels=c("E","F","G","H"),label_y=1.1,label_size=10)

aligned_plots_C <- align_plots((C_barplot_dose_1+theme(legend.position="none",axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7))),
                             (C_barplot_dose_2+theme(legend.position="none",axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7))),
                             (C_barplot_dose_3+theme(legend.position="none",axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7))),
                             (C_barplot_dose_4+theme(legend.position="none",axis.text.x = element_text(size=7),axis.text.y=element_text(size=7),axis.title.x = element_blank(),axis.title.y = element_text(size=7))),
                             align="v",axis="lr")

vax_plots_C <- plot_grid(aligned_plots_C[[1]],aligned_plots_C[[2]],aligned_plots_C[[3]],aligned_plots_C[[4]],
                       ncol=1, rel_heights = c(1,1,1,1),labels=c("I","J","K","L"),label_y = 1.1,label_size=10)

vax_plots <- plot_grid(vax_plots_BNT,vax_plots_C,nrow=1)

vax_plots_with_legend <- plot_grid(vax_plots,legend,ncol=1,rel_heights = c(1,0.1))
vax_plots_with_legend

figure_1 <- plot_grid(incidence_variant_plot,vax_plots_with_legend,ncol=1,rel_heights=c(1,1.1))


ggsave(file = "Fig_1.pdf", figure_1,width = 180, height = 190,units = "mm")

