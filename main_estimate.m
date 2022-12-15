%%% Purpose: To estimate IAR (unvaccinated), VEs, sensitivities and
%%% specificies, and to select posterior samples for IAR, population
%%% immunity and ascertainment rate calculations

%%% Note 1: Vaccination history and viral load from sewage surveillance
%%% data could not be shared pursuant to confidentiality undertakings to the
%%% Department of Health and the Environmental Protection Department, both
%%% of the Government of the HKSAR.

clear all;
close all;
rng(1)

%%%%%%%%%%%%%%%%%% Assume 7 day delay to VE taking effect (primary analysis) %%%%%%%%%%%%%%%%%%%%%%
%%%%% Read data
dir = 'inference_inputs/';
no_vax = readtable([dir 'no_vax_input_data.csv']);  % data file corresponding to unvaccinated individuals            
B = readtable([dir 'b_input_data.csv']);    % data files corresponding to individuals with one-four homologous doses of BNT (BNT162b2)
BB = readtable([dir 'bb_input_data.csv']);                   
BBB = readtable([dir 'bbb_input_data.csv']);              
BBBB = readtable([dir 'bbbb_input_data.csv']);

S = readtable([dir 's_input_data.csv']);    % data files corresponding to individuals with one-four homologous doses of CoronaVac (Sinovac)
SS = readtable([dir 'ss_input_data.csv']);
SSS = readtable([dir 'sss_input_data.csv']);
SSSS = readtable([dir 'ssss_input_data.csv']);

viral_load = readtable([dir 'proxy.csv']);                    
sens_unvax = readtable([dir 'sens_novax_data.csv']);          
sens_BNT = readtable([dir 'sens_BNT_data.csv']);              
sens_C = readtable([dir 'sens_Sinovac_data.csv']);

spec_unvax = readtable([dir 'spec_novax_data.csv']);
spec_BNT = readtable([dir 'spec_BNT_data.csv']);
spec_C = readtable([dir 'spec_C_data.csv']);

%%%%%% Regularise proxies
%%%%%% Start date is set to Jan 1, 2022
proxy = viral_load.two_day_viral_load';
proxy = proxy/sum(proxy);
t_start = 366;
t_end = numel(proxy);
proxy = proxy(t_start:t_end);

%%%%%%  Process vaccination and serum donation date data
num_samples = [size(no_vax,1) size(B,1) size(BB,1) size(BBB,1) size(BBBB,1) size(S,1) size(SS,1) size(SSS,1) size(SSSS,1)];

% Increase vaccination date by the number of days delayed
vax_delay = 7;

% Vax_idx: Group 1 is unvacinated; Group 2 is BNT vaccinated; Group 3 is
% CoronaVac vaccinated
vax_idx = 1;
sero{vax_idx} = no_vax.sero; % seropositivity            
age{vax_idx} = no_vax.age; % participant age
gender{vax_idx} = categorical(no_vax.gender); % participant gender (not used)
donation_day{vax_idx} = no_vax.donation_day-t_start+1; % day when participant donated serum

vax_idx = 2;
sero{vax_idx} = [B.sero; BB.sero; BBB.sero; BBBB.sero];     
age{vax_idx} = [B.age; BB.age; BBB.age; BBBB.age];
gender{vax_idx} = [categorical(B.gender); categorical(BB.gender); categorical(BBB.gender); categorical(BBBB.gender)];
donation_day{vax_idx} = [B.donation_day; BB.donation_day; BBB.donation_day; BBBB.donation_day]-t_start+1;

% Set untaken doses (by serum donation day) to day 9999
B.dose_2_day = repmat(9999,size(B,1),1);
B.dose_3_day = repmat(9999,size(B,1),1);
B.dose_4_day = repmat(9999,size(B,1),1);

BB.dose_3_day = repmat(9999,size(BB,1),1);
BB.dose_4_day = repmat(9999,size(BB,1),1);

BBB.dose_4_day = repmat(9999,size(BBB,1),1);

% Increase vaccination date by the number of days delayed until VE takes effect
vax_delay = 7;

dose_day{vax_idx,1} = [B.dose_1_day; BB.dose_1_day;BBB.dose_1_day;BBBB.dose_1_day]+vax_delay;
dose_day{vax_idx,2} = [B.dose_2_day; BB.dose_2_day;BBB.dose_2_day;BBBB.dose_2_day]+vax_delay;
dose_day{vax_idx,3} = [B.dose_3_day; BB.dose_3_day;BBB.dose_3_day;BBBB.dose_3_day]+vax_delay;
dose_day{vax_idx,4} = [B.dose_4_day; BB.dose_4_day;BBB.dose_4_day;BBBB.dose_4_day]+vax_delay;

vax_idx = 3;
sero{vax_idx} = [S.sero; SS.sero; SSS.sero; SSSS.sero];      
age{vax_idx} = [S.age; SS.age; SSS.age; SSSS.age];
gender{vax_idx} = [categorical(S.gender); categorical(SS.gender); categorical(SSS.gender); categorical(SSSS.gender)];
donation_day{vax_idx} = [S.donation_day; SS.donation_day; SSS.donation_day; SSSS.donation_day]-t_start+1;     

S.dose_2_day = repmat(9999,size(S,1),1);
S.dose_3_day = repmat(9999,size(S,1),1);
S.dose_4_day = repmat(9999,size(S,1),1);

SS.dose_3_day = repmat(9999,size(SS,1),1);
SS.dose_4_day = repmat(9999,size(SS,1),1);

SSS.dose_4_day = repmat(9999,size(SSS,1),1);

dose_day{vax_idx,1} = [S.dose_1_day; SS.dose_1_day;SSS.dose_1_day;SSSS.dose_1_day]+vax_delay;
dose_day{vax_idx,2} = [S.dose_2_day; SS.dose_2_day;SSS.dose_2_day;SSSS.dose_2_day]+vax_delay;
dose_day{vax_idx,3} = [S.dose_3_day; SS.dose_3_day;SSS.dose_3_day;SSSS.dose_3_day]+vax_delay;
dose_day{vax_idx,4} = [S.dose_4_day; SS.dose_4_day;SSS.dose_4_day;SSSS.dose_4_day]+vax_delay;

% Set age PCHIP (piecewise cubic hermite interpolating polynomial) knots for age effect estimation
age_knots = [10 18 35 50 65];
age_knots_7delay = [10 18 35 50 65];


% Obtain and process positive and negative controls data for joint estimation of assay sensitivity
% and specificity 
n_sens_unvax = size(sens_unvax,1);
y_sens_unvax = sum(sens_unvax.sero,1);

n_sens_BNT = size(sens_BNT,1);
y_sens_BNT = sum(sens_BNT.sero,1);

n_sens_Sinovac = size(sens_C,1);
y_sens_Sinovac = sum(sens_C.sero,1);

n_spec_unvax = size(spec_unvax,1);
y_spec_unvax = n_spec_unvax-sum(spec_unvax.sero,1);

n_spec_BNT = size(spec_BNT,1);
y_spec_BNT = n_spec_BNT-sum(spec_BNT.sero,1);

n_spec_Sinovac = size(spec_C,1);
y_spec_Sinovac = n_spec_Sinovac-sum(spec_C.sero,1);


tv = t_start:t_end;
%%%%%% Calculate matrix for tabulating days since the most recent dose for each study
%%%%%% participant at each time t, given delay to VE taking effect


day_since_dose=cell(3,4);
for ii = 2:3
    for jj = 1:4
        day_since_dose{ii,jj} = repmat(tv,numel(dose_day{ii,jj}),1)-repmat(dose_day{ii,jj},1,numel(tv))+1;
    end
end

%%%%%% Calculate matrix containing indicator variables indicating the latest vaccine type and
%%%%%% dose for each study participant at each time t, given delay to
%%%%%% VE taking effect

VE_indicator=cell(3,4);
for ii = 2:3
    for jj = 1:4
            VE_indicator{ii,jj} = day_since_dose{ii,1}>=0;
            for kk = 2:4
                VE_indicator{ii,jj} = (VE_indicator {ii,jj} & day_since_dose{ii,kk}<0)*(jj<kk)+(VE_indicator {ii,jj} & day_since_dose{ii,kk}>=0)*(jj>=kk); 
            end
    end
    proxy_vax{ii} = repmat(proxy,size(day_since_dose{ii,1},1),1);
end


proxy_vax{1} = repmat(proxy,size(sero{1},1),1);


%%% Gather data for inference
data.sero = sero;
data.donation_day = donation_day;
data.proxy_vax = proxy_vax;
data.day_since_dose = day_since_dose;
data.sens = [y_sens_unvax n_sens_unvax; y_sens_BNT n_sens_BNT; y_sens_Sinovac n_sens_Sinovac];
data.spec = [y_spec_unvax n_spec_unvax; y_spec_BNT n_spec_BNT; y_spec_Sinovac n_spec_Sinovac];
data.age = age;
data.age_knots = age_knots;
data.gender = gender;
data.VE_indicator = VE_indicator;

%%%%%%%%%% Perform estimation (see Methods)

chain_7delay.data = data;
chain_7delay.numStepsPerParameter = 2000;
chain_7delay.minNumStepsBeforeRestart = 1000;
chain_7delay.numStepsBetweenDisplay = 1000;
chain_7delay.probAcceptRange = [0.2 0.7];
chain_7delay.UB = [1 repmat(1,1,6) repmat(1,1,2) repmat(1,1,6) repmat(10, 1, numel(age_knots)-1)];      % limit UB for waning rates to 0.2
chain_7delay.numParameters = numel(chain_7delay.UB);
chain_7delay.LB = zeros(1,chain_7delay.numParameters);
chain_7delay.startingPoint = [0.7 repmat(0.5,1,6) repmat(0.01,1,2) repmat(0.9,1,6) ones(1,4)];
chain_7delay.totalNumSteps = chain_7delay.numStepsPerParameter*chain_7delay.numParameters;
chain_7delay.stepSTD = repmat(1, 1, chain_7delay.numParameters);
chain_7delay.logSpace = zeros(1, chain_7delay.numParameters);
chain_7delay.dir = 'Test';

output = logLikelihood(chain_7delay.startingPoint,data);


chain_7delay = MCMC_Gibbs(@logLikelihood, @logPriorDist, chain_7delay);
disp(' ');

posterior_7delay = chain_7delay.record(chain_7delay.totalNumSteps/3:2:end,:);
save('posterior_7delay_age.mat','posterior_7delay');

% Calculate initial VE for 2 or 3 doses of BNT or CoronaVac
posterior_7delay(:,2) = prod(posterior_7delay(:,2:4),2);
posterior_7delay(:,3) = prod(posterior_7delay(:,3:4),2); 

posterior_7delay(:,5) = prod(posterior_7delay(:,5:7),2); 
posterior_7delay(:,6) = prod(posterior_7delay(:,6:7),2); 

% Display posterior medians and 95% credible intervals
var_text = {'IAR_unvax','VE_BNT2','VE_BNT3', 'VE_BNT4','VE_C2','VE_C3','VE_C4','waning_BNT','waning_C','sens_unvax','sens_BNT','sens_C','spec_unvax','spec_BNT','spec_C','10','18','50','65'};
for var_idx = 1:numel(var_text)
    disp(var_text{var_idx});
    disp(prctile(posterior_7delay(:,var_idx), [50 2.5 97.5]));
end

% Reproduce age effect factors using PCHIP
age_range_7delay = min(age_knots_7delay):max(age_knots_7delay);

age_knots_factors_7delay = [prctile(posterior_7delay(:,16:17),[2.5;50;97.5]) ones(3,1) prctile(posterior_7delay(:,18:18+numel(age_knots_7delay)-4),[2.5; 50; 97.5])];

FOI_7delay = 1*pchip(age_knots_7delay,age_knots_factors_7delay,age_range_7delay);

% Calculate VE over time assuming exponential decay (see Methods). VE takes effect on the 7th day.
tv = 0:200;
tv_7delay = 0:(max(tv)-vax_delay);

VE_BNT2_TS_7delay = repmat(posterior_7delay(:,2), 1, numel(tv)).*[zeros(size(posterior_7delay,1),7) exp(-posterior_7delay(:,8)*tv_7delay)];
VE_BNT3_TS_7delay = repmat(posterior_7delay(:,3), 1, numel(tv)).*[zeros(size(posterior_7delay,1),7) exp(-posterior_7delay(:,8)*tv_7delay)];
VE_BNT4_TS_7delay = repmat(posterior_7delay(:,4), 1, numel(tv)).*[zeros(size(posterior_7delay,1),7) exp(-posterior_7delay(:,8)*tv_7delay)];

VE_BNT_prctile_7delay{2} = prctile(VE_BNT2_TS_7delay, [2.5 25 50 75 97.5]);
VE_BNT_prctile_7delay{3} = prctile(VE_BNT3_TS_7delay, [2.5 25 50 75 97.5]);
VE_BNT_prctile_7delay{4} = prctile(VE_BNT4_TS_7delay, [2.5 25 50 75 97.5]);% 

VE_C2_TS_7delay = repmat(posterior_7delay(:,5), 1, numel(tv)).*[zeros(size(posterior_7delay,1),7) exp(-posterior_7delay(:,9)*tv_7delay)];
VE_C3_TS_7delay = repmat(posterior_7delay(:,6), 1, numel(tv)).*[zeros(size(posterior_7delay,1),7) exp(-posterior_7delay(:,9)*tv_7delay)];
VE_C4_TS_7delay = repmat(posterior_7delay(:,7), 1, numel(tv)).*[zeros(size(posterior_7delay,1),7) exp(-posterior_7delay(:,9)*tv_7delay)];

VE_C_prctile_7delay{2} = prctile(VE_C2_TS_7delay, [2.5 25 50 75 97.5]);
VE_C_prctile_7delay{3} = prctile(VE_C3_TS_7delay, [2.5 25 50 75 97.5]);
VE_C_prctile_7delay{4} = prctile(VE_C4_TS_7delay, [2.5 25 50 75 97.5]);% 

save('BNT_C_variables_age.mat');

%%%%%%%%%%%%%%%%%% Assume 14 day delay to VE taking effect (sensitivity analysis) %%%%%%%%%%%%%%%%%%%%%%

% Increase vaccination date by the number of days delayed
vax_delay = 14;

% Notwithstanding the 14 day delay to VE taking effect, the code below follows the same logic as above

vax_idx = 2;
dose_day{vax_idx,1} = [B.dose_1_day; BB.dose_1_day;BBB.dose_1_day;BBBB.dose_1_day]+vax_delay;
dose_day{vax_idx,2} = [B.dose_2_day; BB.dose_2_day;BBB.dose_2_day;BBBB.dose_2_day]+vax_delay;
dose_day{vax_idx,3} = [B.dose_3_day; BB.dose_3_day;BBB.dose_3_day;BBBB.dose_3_day]+vax_delay;
dose_day{vax_idx,4} = [B.dose_4_day; BB.dose_4_day;BBB.dose_4_day;BBBB.dose_4_day]+vax_delay;

vax_idx = 3;
dose_day{vax_idx,1} = [S.dose_1_day; SS.dose_1_day;SSS.dose_1_day;SSSS.dose_1_day]+vax_delay;
dose_day{vax_idx,2} = [S.dose_2_day; SS.dose_2_day;SSS.dose_2_day;SSSS.dose_2_day]+vax_delay;
dose_day{vax_idx,3} = [S.dose_3_day; SS.dose_3_day;SSS.dose_3_day;SSSS.dose_3_day]+vax_delay;
dose_day{vax_idx,4} = [S.dose_4_day; SS.dose_4_day;SSS.dose_4_day;SSSS.dose_4_day]+vax_delay;

age_knots_14delay = [10 18 35 50 65];

n_sens_unvax = size(sens_unvax,1);
y_sens_unvax = sum(sens_unvax.sero,1);

n_sens_BNT = size(sens_BNT,1);
y_sens_BNT = sum(sens_BNT.sero,1);

n_sens_Sinovac = size(sens_C,1);
y_sens_Sinovac = sum(sens_C.sero,1);

n_spec_unvax = size(spec_unvax,1);
y_spec_unvax = n_spec_unvax-sum(spec_unvax.sero,1);

n_spec_BNT = size(spec_BNT,1);
y_spec_BNT = n_spec_BNT-sum(spec_BNT.sero,1);

n_spec_Sinovac = size(spec_C,1);
y_spec_Sinovac = n_spec_Sinovac-sum(spec_C.sero,1);

tv = t_start:t_end;

%%%%%% Recalculate days since last dose matrix given delay

day_since_dose=cell(3,4);
for ii = 2:3
    for jj = 1:4
        day_since_dose{ii,jj} = repmat(tv,numel(dose_day{ii,jj}),1)-repmat(dose_day{ii,jj},1,numel(tv))+1;
    end
end

%%%%%% Recalculate vaccine type and dose indicator variable matrix given delay

VE_indicator=cell(3,4);
for ii = 2:3
    for jj = 1:4
            VE_indicator{ii,jj} = day_since_dose{ii,1}>=0;
            for kk = 2:4
                VE_indicator{ii,jj} = (VE_indicator {ii,jj} & day_since_dose{ii,kk}<0)*(jj<kk)+(VE_indicator {ii,jj} & day_since_dose{ii,kk}>=0)*(jj>=kk); 
            end
    end
    proxy_vax{ii} = repmat(proxy,size(day_since_dose{ii,1},1),1);
end


proxy_vax{1} = repmat(proxy,size(sero{1},1),1);

%%% Gather data for inference
data.sero = sero;
data.donation_day = donation_day;
data.proxy_vax = proxy_vax;
data.day_since_dose = day_since_dose;
data.sens = [y_sens_unvax n_sens_unvax; y_sens_BNT n_sens_BNT; y_sens_Sinovac n_sens_Sinovac];
data.spec = [y_spec_unvax n_spec_unvax; y_spec_BNT n_spec_BNT; y_spec_Sinovac n_spec_Sinovac];
data.age = age;
data.age_knots = age_knots;
data.gender = gender;
data.VE_indicator = VE_indicator;

%%%%%%%%%% Perform inference (see Methods)

chain_14delay.data = data;
chain_14delay.numStepsPerParameter = 2000;
chain_14delay.minNumStepsBeforeRestart = 1000;
chain_14delay.numStepsBetweenDisplay = 1000;
chain_14delay.probAcceptRange = [0.2 0.8];
chain_14delay.UB = [1 repmat(1,1,6) repmat(1,1,2) repmat(1,1,6) repmat(10, 1, numel(age_knots)-1)];      % limit UB for waning rates to 0.2
chain_14delay.numParameters = numel(chain_14delay.UB);
chain_14delay.LB = zeros(1,chain_14delay.numParameters);
chain_14delay.startingPoint = [0.7 repmat(0.5,1,6) repmat(0.01,1,2) repmat(0.9,1,6) ones(1,4)];
chain_14delay.totalNumSteps = chain_14delay.numStepsPerParameter*chain_14delay.numParameters;
chain_14delay.stepSTD = repmat(0.1, 1, chain_14delay.numParameters);
chain_14delay.logSpace = zeros(1, chain_14delay.numParameters);
chain_14delay.dir = 'Test';

output = logLikelihood(chain_14delay.startingPoint,data);


chain_14delay = MCMC_Gibbs(@logLikelihood, @logPriorDist, chain_14delay);
disp(' ');

posterior_14delay = chain_14delay.record(chain_14delay.totalNumSteps/3:2:end,:);
save('posterior_14delay_age.mat','posterior_14delay');

% Calculate initial VE for 2 or 3 doses of BNT or CoronaVac
posterior_14delay(:,2) = prod(posterior_14delay(:,2:4),2);
posterior_14delay(:,3) = prod(posterior_14delay(:,3:4),2); 

posterior_14delay(:,5) = prod(posterior_14delay(:,5:7),2); 
posterior_14delay(:,6) = prod(posterior_14delay(:,6:7),2); 

% Display posterior medians and 95% credible intervals
var_text = {'IAR_unvax','VE_BNT2','VE_BNT3', 'VE_BNT4','VE_C2','VE_C3','VE_C4','waning_BNT','waning_C','sens_unvax','sens_BNT','sens_C','spec_unvax','spec_BNT','spec_C','10','18','50','65'};
for var_idx = 1:numel(var_text)
    disp(var_text{var_idx});
    disp(prctile(posterior_14delay(:,var_idx), [50 2.5 97.5]));
end

% Calculate age effect factors using PCHIP
age_range_14delay = min(age_knots_14delay):max(age_knots_14delay);

age_knots_factors_14delay = [prctile(posterior_14delay(:,16:17),[2.5;50;97.5]) ones(3,1) prctile(posterior_14delay(:,18:18+numel(age_knots_14delay)-4),[2.5; 50; 97.5])];

FOI_14delay = 1*pchip(age_knots_14delay,age_knots_factors_14delay,age_range_14delay);

% Calculate VE over time assuming exponential decay (see Methods). VE takes effect on the 14th day.
tv = 0:200;
tv_14delay = 0:(max(tv)-vax_delay);

VE_BNT2_TS_14delay = repmat(posterior_14delay(:,2), 1, numel(tv)).*[zeros(size(posterior_14delay,1),vax_delay) exp(-posterior_14delay(:,8)*tv_14delay)];
VE_BNT3_TS_14delay = repmat(posterior_14delay(:,3), 1, numel(tv)).*[zeros(size(posterior_14delay,1),vax_delay) exp(-posterior_14delay(:,8)*tv_14delay)];
VE_BNT4_TS_14delay = repmat(posterior_14delay(:,4), 1, numel(tv)).*[zeros(size(posterior_14delay,1),vax_delay) exp(-posterior_14delay(:,8)*tv_14delay)];

VE_BNT_prctile_14delay{2} = prctile(VE_BNT2_TS_14delay, [2.5 25 50 75 97.5]);
VE_BNT_prctile_14delay{3} = prctile(VE_BNT3_TS_14delay, [2.5 25 50 75 97.5]);
VE_BNT_prctile_14delay{4} = prctile(VE_BNT4_TS_14delay, [2.5 25 50 75 97.5]);% 

VE_C2_TS_14delay = repmat(posterior_14delay(:,5), 1, numel(tv)).*[zeros(size(posterior_14delay,1),vax_delay) exp(-posterior_14delay(:,9)*tv_14delay)];
VE_C3_TS_14delay = repmat(posterior_14delay(:,6), 1, numel(tv)).*[zeros(size(posterior_14delay,1),vax_delay) exp(-posterior_14delay(:,9)*tv_14delay)];
VE_C4_TS_14delay = repmat(posterior_14delay(:,7), 1, numel(tv)).*[zeros(size(posterior_14delay,1),vax_delay) exp(-posterior_14delay(:,9)*tv_14delay)];

VE_C_prctile_14delay{2} = prctile(VE_C2_TS_14delay, [2.5 25 50 75 97.5]);
VE_C_prctile_14delay{3} = prctile(VE_C3_TS_14delay, [2.5 25 50 75 97.5]);
VE_C_prctile_14delay{4} = prctile(VE_C4_TS_14delay, [2.5 25 50 75 97.5]);% 

save('BNT_C_variables_age.mat');

%%%%%%%%%%%%%%%%%% Assume 21 day delay to VE taking effect (sensitivity analysis) %%%%%%%%%%%%%%%%%%%%%%

% Increase vaccination date by the number of days delayed
vax_delay = 21;

vax_idx = 2;
dose_day{vax_idx,1} = [B.dose_1_day; BB.dose_1_day;BBB.dose_1_day;BBBB.dose_1_day]+vax_delay;
dose_day{vax_idx,2} = [B.dose_2_day; BB.dose_2_day;BBB.dose_2_day;BBBB.dose_2_day]+vax_delay;
dose_day{vax_idx,3} = [B.dose_3_day; BB.dose_3_day;BBB.dose_3_day;BBBB.dose_3_day]+vax_delay;
dose_day{vax_idx,4} = [B.dose_4_day; BB.dose_4_day;BBB.dose_4_day;BBBB.dose_4_day]+vax_delay;

vax_idx = 3;
dose_day{vax_idx,1} = [S.dose_1_day; SS.dose_1_day;SSS.dose_1_day;SSSS.dose_1_day]+vax_delay;
dose_day{vax_idx,2} = [S.dose_2_day; SS.dose_2_day;SSS.dose_2_day;SSSS.dose_2_day]+vax_delay;
dose_day{vax_idx,3} = [S.dose_3_day; SS.dose_3_day;SSS.dose_3_day;SSSS.dose_3_day]+vax_delay;
dose_day{vax_idx,4} = [S.dose_4_day; SS.dose_4_day;SSS.dose_4_day;SSSS.dose_4_day]+vax_delay;

age_knots_21delay = [10 18 35 50 65];

n_sens_unvax = size(sens_unvax,1);
y_sens_unvax = sum(sens_unvax.sero,1);

n_sens_BNT = size(sens_BNT,1);
y_sens_BNT = sum(sens_BNT.sero,1);

n_sens_Sinovac = size(sens_C,1);
y_sens_Sinovac = sum(sens_C.sero,1);

n_spec_unvax = size(spec_unvax,1);
y_spec_unvax = n_spec_unvax-sum(spec_unvax.sero,1);

n_spec_BNT = size(spec_BNT,1);
y_spec_BNT = n_spec_BNT-sum(spec_BNT.sero,1);

n_spec_Sinovac = size(spec_C,1);
y_spec_Sinovac = n_spec_Sinovac-sum(spec_C.sero,1);

tv = t_start:t_end;

%%%%%% Recalculate days since last dose matrix given delay

day_since_dose=cell(3,4);
for ii = 2:3
    for jj = 1:4
        day_since_dose{ii,jj} = repmat(tv,numel(dose_day{ii,jj}),1)-repmat(dose_day{ii,jj},1,numel(tv))+1;
    end
end

%%%%%% Recalculate vaccine type and dose indicator variable matrix given delay

VE_indicator=cell(3,4);
for ii = 2:3
    for jj = 1:4
            VE_indicator{ii,jj} = day_since_dose{ii,1}>=0;
            for kk = 2:4
                VE_indicator{ii,jj} = (VE_indicator {ii,jj} & day_since_dose{ii,kk}<0)*(jj<kk)+(VE_indicator {ii,jj} & day_since_dose{ii,kk}>=0)*(jj>=kk); 
            end
    end
    proxy_vax{ii} = repmat(proxy,size(day_since_dose{ii,1},1),1);
end


proxy_vax{1} = repmat(proxy,size(sero{1},1),1);

%%% Gather data for inference
data.sero = sero;
data.donation_day = donation_day;
data.proxy_vax = proxy_vax;
data.day_since_dose = day_since_dose;
data.sens = [y_sens_unvax n_sens_unvax; y_sens_BNT n_sens_BNT; y_sens_Sinovac n_sens_Sinovac];
data.spec = [y_spec_unvax n_spec_unvax; y_spec_BNT n_spec_BNT; y_spec_Sinovac n_spec_Sinovac];
data.age = age;
data.age_knots = age_knots;
data.gender = gender;
data.VE_indicator = VE_indicator;

%%%%%%%%%% Perform inference (see Methods)

chain_21delay.data = data;
chain_21delay.numStepsPerParameter = 2000;
chain_21delay.minNumStepsBeforeRestart = 1000;
chain_21delay.numStepsBetweenDisplay = 1000;
chain_21delay.probAcceptRange = [0.2 0.8];
chain_21delay.UB = [1 repmat(1,1,6) repmat(1,1,2) repmat(1,1,6) repmat(10, 1, numel(age_knots)-1)];      % limit UB for waning rates to 0.2
chain_21delay.numParameters = numel(chain_21delay.UB);
chain_21delay.LB = zeros(1,chain_21delay.numParameters);
chain_21delay.startingPoint = [0.7 repmat(0.5,1,6) repmat(0.01,1,2) repmat(0.9,1,6) ones(1,4)];
chain_21delay.totalNumSteps = chain_21delay.numStepsPerParameter*chain_21delay.numParameters;
chain_21delay.stepSTD = repmat(0.1, 1, chain_21delay.numParameters);
chain_21delay.logSpace = zeros(1, chain_21delay.numParameters);
chain_21delay.dir = 'Test';

output = logLikelihood(chain_21delay.startingPoint,data);


chain_21delay = MCMC_Gibbs(@logLikelihood, @logPriorDist, chain_21delay);
disp(' ');

posterior_21delay = chain_21delay.record(chain_21delay.totalNumSteps/3:2:end,:);
save('posterior_21delay_age.mat','posterior_21delay');

% Calculate initial VE for 2 or 3 doses of BNT or CoronaVac
posterior_21delay(:,2) = prod(posterior_21delay(:,2:4),2);
posterior_21delay(:,3) = prod(posterior_21delay(:,3:4),2); 

posterior_21delay(:,5) = prod(posterior_21delay(:,5:7),2); 
posterior_21delay(:,6) = prod(posterior_21delay(:,6:7),2); 

% Display posterior medians and 95% credible intervals
var_text = {'IAR_unvax','VE_BNT2','VE_BNT3', 'VE_BNT4','VE_C2','VE_C3','VE_C4','waning_BNT','waning_C','sens_unvax','sens_BNT','sens_C','spec_unvax','spec_BNT','spec_C','10','18','50','65'};
for var_idx = 1:numel(var_text)
    disp(var_text{var_idx});
    disp(prctile(posterior_21delay(:,var_idx), [50 2.5 97.5]));
end

% Reproduce age effect factors using PCHIP
age_range_21delay = min(age_knots_21delay):max(age_knots_21delay);

age_knots_factors_21delay = [prctile(posterior_21delay(:,16:17),[2.5;50;97.5]) ones(3,1) prctile(posterior_21delay(:,18:18+numel(age_knots_21delay)-4),[2.5; 50; 97.5])];

FOI_21delay = 1*pchip(age_knots_21delay,age_knots_factors_21delay,age_range_21delay);

% Calculate VE over time assuming exponential decay (see Methods). VE takes effect on the 21st day.
tv = 0:200;
tv_21delay = 0:(max(tv)-vax_delay);

VE_BNT2_TS_21delay = repmat(posterior_21delay(:,2), 1, numel(tv)).*[zeros(size(posterior_21delay,1),vax_delay) exp(-posterior_21delay(:,8)*tv_21delay)];
VE_BNT3_TS_21delay = repmat(posterior_21delay(:,3), 1, numel(tv)).*[zeros(size(posterior_21delay,1),vax_delay) exp(-posterior_21delay(:,8)*tv_21delay)];
VE_BNT4_TS_21delay = repmat(posterior_21delay(:,4), 1, numel(tv)).*[zeros(size(posterior_21delay,1),vax_delay) exp(-posterior_21delay(:,8)*tv_21delay)];

VE_BNT_prctile_21delay{2} = prctile(VE_BNT2_TS_21delay, [2.5 25 50 75 97.5]);
VE_BNT_prctile_21delay{3} = prctile(VE_BNT3_TS_21delay, [2.5 25 50 75 97.5]);
VE_BNT_prctile_21delay{4} = prctile(VE_BNT4_TS_21delay, [2.5 25 50 75 97.5]);% 

VE_C2_TS_21delay = repmat(posterior_21delay(:,5), 1, numel(tv)).*[zeros(size(posterior_21delay,1),vax_delay) exp(-posterior_21delay(:,9)*tv_21delay)];
VE_C3_TS_21delay = repmat(posterior_21delay(:,6), 1, numel(tv)).*[zeros(size(posterior_21delay,1),vax_delay) exp(-posterior_21delay(:,9)*tv_21delay)];
VE_C4_TS_21delay = repmat(posterior_21delay(:,7), 1, numel(tv)).*[zeros(size(posterior_21delay,1),vax_delay) exp(-posterior_21delay(:,9)*tv_21delay)];

VE_C_prctile_21delay{2} = prctile(VE_C2_TS_21delay, [2.5 25 50 75 97.5]);
VE_C_prctile_21delay{3} = prctile(VE_C3_TS_21delay, [2.5 25 50 75 97.5]);
VE_C_prctile_21delay{4} = prctile(VE_C4_TS_21delay, [2.5 25 50 75 97.5]);% 

save('BNT_C_variables_age.mat');


%%% Plot Extended Data Figure 1
figure(1);
clf;
hold on;
dat_7delay = FOI_7delay;
dat_14delay = FOI_14delay;
dat_21delay = FOI_21delay;

h_7delay = fill([age_range fliplr(age_range)], [dat_7delay(1,:) fliplr(dat_7delay(end,:))], [0.902 0.294 0.208]);
transparency = 0.2;
set(h_7delay, "FaceAlpha", transparency, "lineStyle", "none");

h_14delay = fill([age_range fliplr(age_range)], [dat_14delay(1,:) fliplr(dat_14delay(end,:))], [0.302 0.733 1]);
transparency = 0.2;
set(h_14delay, "FaceAlpha", transparency, "lineStyle", "none");

h_21delay = fill([age_range fliplr(age_range)], [dat_21delay(1,:) fliplr(dat_21delay(end,:))], "g");
transparency = 0.2;
set(h_21delay, "FaceAlpha", transparency, "lineStyle", "none");


plot(age_range, dat_7delay(round(size(dat_7delay,1)/2),:),"Color",[0.902 0.294 0.208]);
plot(age_range, dat_14delay(round(size(dat_14delay,1)/2),:),"Color",[0.302 0.733 1]);
plot(age_range, dat_21delay(round(size(dat_21delay,1)/2),:),"Color","g");


axis tight;
xlabel("Age");
ylabel("Relative FOI vs 35 year old subject");
legend("7 day delay","14 day delay","21 day delay");

exportgraphics(gcf,'extended_data_fig_1.tiff','Resolution',300);

save("Extended_data_Figure_1.mat",'age_range','dat_7delay','dat_14delay','dat_21delay');


%%% Randomly select 300 posterior samples for IAR and population immunity
%%% calculation
sampled_posteriors_7delay_dir = 'IAR_calculation_7delay/';
sampled_posteriors_14delay_dir = 'IAR_calculation_14delay/';
sampled_posteriors_21delay_dir = 'IAR_calculation_21delay/';

num_sampled_posteriors_7delay = randsample([1:size(posterior_7delay,1)],300);
sampled_posteriors_7delay = posterior_7delay(num_sampled_posteriors_7delay,:);

save([sampled_posteriors_7delay_dir 'sampled_posteriors_7delay_age.mat'],'sampled_posteriors_7delay')

num_sampled_posteriors_14delay = randsample([1:size(posterior_14delay,1)],300);
sampled_posteriors_14delay = posterior_14delay(num_sampled_posteriors_14delay,:);

save([sampled_posteriors_14delay_dir 'sampled_posteriors_14delay_age.mat'],'sampled_posteriors_14delay')

num_sampled_posteriors_21delay = randsample([1:size(posterior_21delay,1)],300);
sampled_posteriors_21delay = posterior_21delay(num_sampled_posteriors_21delay,:);

save([sampled_posteriors_21delay_dir 'sampled_posteriors_21delay_age.mat'],'sampled_posteriors_21delay')



%%% Output VE over time for creating Figure 2
output_dir = 'Figure_2/VE/';
writematrix([tv; VE_BNT_prctile_7delay{2}]',[output_dir 'VE_BB_prctile_7delay_age.csv'])
writematrix([tv; VE_BNT_prctile_7delay{3}]',[output_dir 'VE_BBB_prctile_7delay_age.csv'])
writematrix([tv; VE_BNT_prctile_7delay{4}]',[output_dir 'VE_BBBB_prctile_7delay_age.csv'])

writematrix([tv; VE_C_prctile_7delay{2}]',[output_dir 'VE_CC_prctile_7delay_age.csv'])
writematrix([tv; VE_C_prctile_7delay{3}]',[output_dir 'VE_CCC_prctile_7delay_age.csv'])
writematrix([tv; VE_C_prctile_7delay{4}]',[output_dir 'VE_CCCC_prctile_7delay_age.csv'])

writematrix([tv; VE_BNT_prctile_14delay{2}]',[output_dir 'VE_BB_prctile_14delay_age.csv'])
writematrix([tv; VE_BNT_prctile_14delay{3}]',[output_dir 'VE_BBB_prctile_14delay_age.csv'])
writematrix([tv; VE_BNT_prctile_14delay{4}]',[output_dir 'VE_BBBB_prctile_14delay_age.csv'])

writematrix([tv; VE_C_prctile_14delay{2}]',[output_dir 'VE_CC_prctile_14delay_age.csv'])
writematrix([tv; VE_C_prctile_14delay{3}]',[output_dir 'VE_CCC_prctile_14delay_age.csv'])
writematrix([tv; VE_C_prctile_14delay{4}]',[output_dir 'VE_CCCC_prctile_14delay_age.csv'])

writematrix([tv; VE_BNT_prctile_21delay{2}]',[output_dir 'VE_BB_prctile_21delay_age.csv'])
writematrix([tv; VE_BNT_prctile_21delay{3}]',[output_dir 'VE_BBB_prctile_21delay_age.csv'])
writematrix([tv; VE_BNT_prctile_21delay{4}]',[output_dir 'VE_BBBB_prctile_21delay_age.csv'])

writematrix([tv; VE_C_prctile_21delay{2}]',[output_dir 'VE_CC_prctile_21delay_age.csv'])
writematrix([tv; VE_C_prctile_21delay{3}]',[output_dir 'VE_CCC_prctile_21delay_age.csv'])
writematrix([tv; VE_C_prctile_21delay{4}]',[output_dir 'VE_CCCC_prctile_21delay_age.csv'])

