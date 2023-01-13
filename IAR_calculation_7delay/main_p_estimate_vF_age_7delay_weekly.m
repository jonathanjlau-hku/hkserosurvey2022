clear all;
close all;

% Purpose: Calculate IAR and population immunity assuming VE takes
% effect with a 7 day delay and given no waning of infection-acquired
% immunity, 15% waning over 365 days, or 25% waning over 100 days.

% Note 1: Vaccination records and viral load data from sewage surveillance
% cannot be shared due to confidentiality Undertakings to the Department of Health and the
% Environmental Protection Department, both of the Government of the HKSAR,
% respectively.

c = parcluster('local');
c.NumWorkers = 32;
p = c.parpool(32);

clear("7_delay_age_weekly.mat")

%%%%%% Read data
IAR_records = readmatrix("all_records.csv",'NumHeaderLines',1);            
totalPopu = size(IAR_records,1);
VE_records = zeros(size(IAR_records,1),151);
VE_records(:,1) = IAR_records(:,2);
viral_load = readtable("proxy.csv");                   
load("sampled_posteriors_7delay_age.mat");
sampled_posteriors = sampled_posteriors_7delay;

%%%%%% Regularise proxy against start date of January 1, 2022
t_start = 366;
proxy = viral_load.two_day_viral_load';
proxy = proxy/sum(proxy);
t_end = min(numel(proxy),577);
proxy = proxy(t_start:t_end);
cum_proxy = cumsum(proxy);
proxyLength = numel(proxy);

%%%%%% Collect age and vaccination history for each of the 7 million individuals in Hong Kong. 
%%%%%% Set age maximum at 85.
tv = t_start:t_end;
tv_day = tv-t_start+1;
flip_tv_day = flip(tv_day);

tv_week = [t_start:7:t_end]-t_start+1;
tv_week = [tv_week 212];

proxyWeekLength = numel(tv_week);
vax_delay = 7;

age = min(2022-IAR_records(:,2),85);
age_knots = [0 10 18 35 50 65 100];
type_dose_1 = IAR_records(:,3);
day_dose_1 = IAR_records(:,4)+vax_delay-1;
type_dose_2 = IAR_records(:,5);
day_dose_2 = IAR_records(:,6)+vax_delay-1;
type_dose_3 = IAR_records(:,7);
day_dose_3 = IAR_records(:,8)+vax_delay-1;
type_dose_4 = IAR_records(:,9);
day_dose_4 = IAR_records(:,10)+vax_delay-1;
no_dose = IAR_records(:,11);

individuals = numel(type_dose_1);

%%% Infection-induced immunity waning rates based on Barnard(2022), Malato
%%% (2022) and Altarawneh (2022)
waning_barnard = -log(0.85)/365; 
waning_malato = -log(0.75)/100;

%%% Prepare matfile (v7.3)
save ("7_delay_age_weekly.mat","proxy","-v7.3","-nocompression")
m = matfile("7_delay_age_weekly - Copy.mat","writable",true);

m.IAR_records = IAR_records;

p_infection = zeros(totalPopu,proxyWeekLength);
immunity = zeros(totalPopu,proxyWeekLength);
immunity_barnard = zeros(totalPopu,proxyWeekLength);
immunity_malato = zeros(totalPopu,proxyWeekLength);

% Calculate p(infection) and immunity (see Methods) for each individual
% across each of the 300 posterior samples

for ii = [2:-1:1 3:300]
    ii
    IAR_unvax=sampled_posteriors(ii,1);
    VE_BNT1=0;
    VE_BNT2=sampled_posteriors(ii,2);
    VE_BNT3=sampled_posteriors(ii,3);
    VE_BNT4=sampled_posteriors(ii,4);
    VE_C1=0;
    VE_C2=sampled_posteriors(ii,5);
    VE_C3=sampled_posteriors(ii,6);
    VE_C4=sampled_posteriors(ii,7);
    waning_BNT=sampled_posteriors(ii,8);
    waning_C=sampled_posteriors(ii,9);
    age_knots_factors=[sampled_posteriors(ii,16) sampled_posteriors(ii,16:17) 1 sampled_posteriors(ii,18:end) sampled_posteriors(ii,end)];
    age_scaling_factor = pchip(age_knots,age_knots_factors,0:100);
    
    % loop through each individual to calculate prop of infection and
    % immunity
        
    jj_infect = zeros(1,proxyWeekLength);
    jj_immunity = zeros(1,proxyWeekLength);
    jj_immunity_barnard = zeros(1,proxyWeekLength);
    jj_immunity_malato = zeros(1,proxyWeekLength);
    
    ii
    parfor jj = 1:individuals
        
         [jj_infect,jj_immunity,jj_immunity_barnard,jj_immunity_malato] = infection_calc_vF_weekly(IAR_unvax,VE_BNT1,VE_BNT2,VE_BNT3,VE_BNT4,VE_C1,VE_C2,VE_C3,VE_C4,...
            waning_BNT, waning_C,waning_barnard,waning_malato,age_scaling_factor(age(jj)+1),proxy,tv,tv_week,flip_tv_day,...
            type_dose_1(jj),day_dose_1(jj),type_dose_2(jj),day_dose_2(jj),type_dose_3(jj),day_dose_3(jj),...
            type_dose_4(jj),day_dose_4(jj),no_dose(jj));

        p_infection(jj,:) = single(jj_infect);
        immunity(jj,:) = single(jj_immunity);
        immunity_barnard(jj,:) = single(jj_immunity_barnard);
        immunity_malato(jj,:) = single(jj_immunity_malato);
        
    end
    ii
    m.p_infection(1:totalPopu,1:proxyWeekLength,ii) = p_infection;
    m.immunity(1:totalPopu,1:proxyWeekLength,ii) = immunity;
    m.immunity_barnard(1:totalPopu,1:proxyWeekLength,ii) = immunity_barnard;
    m.immunity_malato(1:totalPopu,1:proxyWeekLength,ii) = immunity_malato;


end

age_thresholds = [-1 11 19 29 39 49 59 100];
age_cohorts = ["0-11" "12-19" "20-29" "30-39" "40-49" "50-59" "60+"];

p_infection_weekly_cohort = zeros(numel(age_cohorts),6,proxyWeekLength);
immunity_weekly_cohort = zeros(numel(age_cohorts),6,proxyWeekLength);
immunity_barnard_weekly_cohort = zeros(numel(age_cohorts),6,proxyWeekLength);
immunity_malato_weekly_cohort = zeros(numel(age_cohorts),6,proxyWeekLength);


p_infection_weekly_all = zeros(1,6,proxyWeekLength);
immunity_weekly_all = zeros(1,6,proxyWeekLength);
immunity_barnard_weekly_all = zeros(1,6,proxyWeekLength);
immunity_malato_weekly_all = zeros(1,6,proxyWeekLength);

no_scenarios = 300;

% Calculate weekly IAR and population immunity under each waning assumption
% (median and 95% CI) over the fifth wave from Jan 1, 2022 to July 31, 2022


parfor ii = 1:proxyWeekLength
    ii
    temp_p_infection_cohort = zeros(numel(age_cohorts),no_scenarios);
    temp_immunity_cohort = zeros(numel(age_cohorts),no_scenarios);
    temp_immunity_barnard_cohort = zeros(numel(age_cohorts),no_scenarios);
    temp_immunity_malato_cohort = zeros(numel(age_cohorts),no_scenarios);

    
    temp_p_infection_all = zeros(1,no_scenarios);
    temp_immunity_all = zeros(1,no_scenarios);
    temp_immunity_barnard_all = zeros(1,no_scenarios);
    temp_immunity_malato_all = zeros(1,no_scenarios);
    
    
    for jj = 1:no_scenarios
        jj
        p_infection_iter = m.p_infection(1:individuals,ii,jj);
        immunity_iter = m.immunity(1:individuals,ii,jj);
        immunity_barnard_iter = m.immunity_barnard(1:individuals,ii,jj);
        immunity_malato_iter = m.immunity_malato(1:individuals,ii,jj);

                
        for kk = 1:numel(age_thresholds)-1
            temp_p_infection_cohort(kk,jj) = mean(p_infection_iter(age> age_thresholds(kk) & age<=age_thresholds(kk+1)));
            temp_immunity_cohort(kk,jj) = mean(immunity_iter(age> age_thresholds(kk) & age<= age_thresholds(kk+1)));
            temp_immunity_barnard_cohort(kk,jj) = mean(immunity_barnard_iter(age> age_thresholds(kk) & age<= age_thresholds(kk+1)));
            temp_immunity_malato_cohort(kk,jj) = mean(immunity_malato_iter(age> age_thresholds(kk) & age<= age_thresholds(kk+1)));
                        
        end
        
        temp_p_infection_all(jj) = mean(p_infection_iter);
        temp_immunity_all(jj) = mean(immunity_iter);
        temp_immunity_barnard_all(jj) = mean(immunity_barnard_iter);
        temp_immunity_malato_all(jj) = mean(immunity_malato_iter);
                
    end
        
    p_infection_weekly_cohort(:,:,ii) = [mean(temp_p_infection_cohort,2) prctile(temp_p_infection_cohort,[2.5 25 50 75 97.5],2)];
    immunity_weekly_cohort(:,:,ii) = [mean(temp_immunity_cohort,2) prctile(temp_immunity_cohort,[2.5 25 50 75 97.5],2)];
    immunity_barnard_weekly_cohort(:,:,ii) = [mean(temp_immunity_barnard_cohort,2) prctile(temp_immunity_barnard_cohort,[2.5 25 50 75 97.5],2)];
    immunity_malato_weekly_cohort(:,:,ii) = [mean(temp_immunity_malato_cohort,2) prctile(temp_immunity_malato_cohort,[2.5 25 50 75 97.5],2)];


    p_infection_weekly_all(1,:,ii) = [mean(temp_p_infection_all) prctile(temp_p_infection_all,[2.5 25 50 75 97.5])];
    immunity_weekly_all(1,:,ii) = [mean(temp_immunity_all) prctile(temp_immunity_all,[2.5 25 50 75 97.5])];
    immunity_barnard_weekly_all(1,:,ii) = [mean(temp_immunity_barnard_all) prctile(temp_immunity_barnard_all,[2.5 25 50 75 97.5])];
    immunity_malato_weekly_all(1,:,ii) = [mean(temp_immunity_malato_all) prctile(temp_immunity_malato_all,[2.5 25 50 75 97.5])];
end


% Summarise results and output
for kk = 1:numel(age_thresholds)-1
    
writematrix([reshape(p_infection_weekly_cohort(kk,:,:),[6 proxyWeekLength])' repmat(age_cohorts(kk),proxyWeekLength,1) repmat("7 day delay",proxyWeekLength,1) repmat("IAR",proxyWeekLength,1)],strcat("p_infection_weekly_7delay_",age_cohorts(kk),".csv"));
writematrix([reshape(immunity_weekly_cohort(kk,:,:),[6 proxyWeekLength])' repmat(age_cohorts(kk),proxyWeekLength,1) repmat("7 day delay",proxyWeekLength,1) repmat("Pop immunity",proxyWeekLength,1)],strcat("immunity_weekly_7delay_",age_cohorts(kk),".csv"));
writematrix([reshape(immunity_barnard_weekly_cohort(kk,:,:),[6 proxyWeekLength])' repmat(age_cohorts(kk),proxyWeekLength,1) repmat("7 day delay",proxyWeekLength,1) repmat("Pop immunity (Scenario 1)",proxyWeekLength,1)],strcat("immunity_weekly_7delay_1_",age_cohorts(kk),".csv"));
writematrix([reshape(immunity_malato_weekly_cohort(kk,:,:),[6 proxyWeekLength])' repmat(age_cohorts(kk),proxyWeekLength,1) repmat("7 day delay",proxyWeekLength,1) repmat("Pop immunity (Scenario 2)",proxyWeekLength,1)],strcat("immunity_weekly_7delay_2_",age_cohorts(kk),".csv"));



end

writematrix([reshape(p_infection_weekly_all,[6 proxyWeekLength])' repmat("all",proxyWeekLength,1) repmat("7 day delay",proxyWeekLength,1) repmat("IAR",proxyWeekLength,1)],"p_infection_weekly_7delay_all.csv")
writematrix([reshape(immunity_weekly_all,[6 proxyWeekLength])' repmat("all",proxyWeekLength,1) repmat("7 day delay",proxyWeekLength,1) repmat("Pop immunity",proxyWeekLength,1)],"immunity_weekly_7delay_all.csv")
writematrix([reshape(immunity_barnard_weekly_all,[6 proxyWeekLength])' repmat("all",proxyWeekLength,1) repmat("7 day delay",proxyWeekLength,1) repmat("Pop immunity (Scenario 1)",proxyWeekLength,1)],"immunity_weekly_7delay_1_all.csv")
writematrix([reshape(immunity_malato_weekly_all,[6 proxyWeekLength])' repmat("all",proxyWeekLength,1) repmat("7 day delay",proxyWeekLength,1) repmat("Pop immunity (Scenario 2)",proxyWeekLength,1)],"immunity_weekly_7delay_2_all.csv")

