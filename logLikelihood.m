function funval = logLikelihood(x, data)

%%% Purpose: Function to calculate log likelihood for the purposes of MCMC inference (see Methods) 

%%%%%% data
sero = data.sero;
donation_day = data.donation_day;
day_since_dose = data.day_since_dose;
sens = data.sens;
spec = data.spec;
age = data.age;
age_knots = data.age_knots;
VE_indicator = data.VE_indicator;

%%%%%% parameters
IAR_unvax = x(1);
VE_BNT1 = 0;
VE_BNT2 = prod(x(2:4));
VE_BNT3 = prod(x(3:4));
VE_BNT4 = x(4);
VE_C1 = 0;
VE_C2 = prod(x(5:7));
VE_C3 = prod(x(6:7));
VE_C4 = x(7);
waning_rate_BNT = x(8);
waning_rate_C = x(9);
p_sens = x(10:12);
p_spec = x(13:15);
age_knots_factors = [x(16) x(17) 1 x(18:18+numel(age_knots)-4)];
total_FOI = -log(1-IAR_unvax);

age_knots_factors = [age_knots_factors(1) age_knots_factors age_knots_factors(end)];
age_knots = [0 age_knots 100];

proxy_unvax = data.proxy_vax{1};
cum_proxy_unvax = zeros(numel(sero{1}),1);

% Caclulate cumulative proxy for unvaccinated (see Methods)
for ii = 1:numel(sero{1})
    cum_proxy_unvax(ii) = sum(proxy_unvax(ii,1:donation_day{1}(ii)-1),2);
end

% Calculate via PCHIP and apply age scaling factors for every study subject
% given their age
for ii =1:numel(age)
     age_scaling_factors{ii} = pchip(age_knots, age_knots_factors,age{ii});
end

% Calculate cumulative proxy for BNT vaccinated individuals (see Methods)
VE_mat_BNT = VE_BNT1*exp(-waning_rate_BNT*day_since_dose{2,1}).*VE_indicator{2,1}+...
    VE_BNT2*exp(-waning_rate_BNT*day_since_dose{2,2}).*VE_indicator{2,2}+...
    VE_BNT3*exp(-waning_rate_BNT*day_since_dose{2,3}).*VE_indicator{2,3}+...
    VE_BNT4*exp(-waning_rate_BNT*day_since_dose{2,4}).*VE_indicator{2,4};


proxy_BNT = data.proxy_vax{2}.*(1-VE_mat_BNT);

cum_proxy_BNT = zeros(numel(sero{2}),1);

for ii = 1:numel(sero{2})
    cum_proxy_BNT(ii) = sum(proxy_BNT(ii,1:donation_day{2}(ii)-1),2);
end

% Calculate cumulative proxy for CoronaVac vaccinated individuals (see
% Methods)
VE_mat_C = VE_C1*exp(-waning_rate_C*day_since_dose{3,1}).*VE_indicator{3,1}+...
    VE_C2*exp(-waning_rate_C*day_since_dose{3,2}).*VE_indicator{3,2}+...
    VE_C3*exp(-waning_rate_C*day_since_dose{3,3}).*VE_indicator{3,3}+...
    VE_C4*exp(-waning_rate_C*day_since_dose{3,4}).*VE_indicator{3,4};



proxy_C = data.proxy_vax{3}.*(1-VE_mat_C);

cum_proxy_C = zeros(numel(sero{3}),1);
for ii = 1:numel(sero{3})
    cum_proxy_C(ii) = sum(proxy_C(ii,1:donation_day{3}(ii)-1),2);
end

% Calculate probability of infection for each individual and calculate input into the likelihood function (see Methods)
for vax_idx=1:numel(sero)
    yy = sero{vax_idx};
    if (vax_idx==1)
        p_infection{vax_idx} = 1-exp(-(age_scaling_factors{vax_idx}.*total_FOI).*cum_proxy_unvax)';
    elseif (vax_idx==2)
        p_infection{vax_idx} = 1-exp(-(age_scaling_factors{vax_idx}.*total_FOI).*cum_proxy_BNT)';
    elseif (vax_idx==3)
        p_infection{vax_idx} = 1-exp(-(age_scaling_factors{vax_idx}.*total_FOI).*cum_proxy_C)';
    end
    p_seropos{vax_idx} = p_infection{vax_idx}*p_sens(vax_idx)+(1-p_infection{vax_idx})*(1-p_spec(vax_idx));
    funval(vax_idx) = loge(p_seropos{vax_idx})*yy+loge(1-p_seropos{vax_idx})*(1-yy);
   
    
end

funval = sum(funval);

% Incorporate test sensitivity and specficity in likelihood function
for ii = 1:numel(p_sens)
    funval = funval+sens(ii,1)*loge(p_sens(ii))+(sens(ii,2)-sens(ii,1))*loge(1-p_sens(ii));
end

for ii = 1:numel(p_spec)
    funval = funval+spec(ii,1)*loge(p_spec(ii))+(spec(ii,2)-spec(ii,1))*loge(1-p_spec(ii));
end