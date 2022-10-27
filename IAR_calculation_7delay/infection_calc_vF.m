function [p_infection_week,immunity_week] = infection_calc_vF(IAR_unvax,VE_BNT1,VE_BNT2,VE_BNT3,VE_BNT4,VE_C1,VE_C2,VE_C3,VE_C4,...
        waning_BNT,waning_C,age_scaling_factor,proxy,tv,tv_week,...
        type_dose_1,day_dose_1,type_dose_2,day_dose_2,type_dose_3,day_dose_3,...
        type_dose_4,day_dose_4,no_dose)
    
    
% Purpose: Calculate p_infection and population immunity for each
% individual

% Calculate cumulative force of infection amongst unvaccinated, days since
% dose and vaccine type and dose indivator variable
total_FOI = -log(1-IAR_unvax);

day_since_dose_1 = tv-repmat(day_dose_1,1,numel(tv));
day_since_dose_2 = tv-repmat(day_dose_2,1,numel(tv));
day_since_dose_3 = tv-repmat(day_dose_3,1,numel(tv));
day_since_dose_4 = tv-repmat(day_dose_4,1,numel(tv));

VE_1_indicator = day_since_dose_1>=0 & day_since_dose_2 <0 & day_since_dose_3 <0 & day_since_dose_4 <0;
VE_2_indicator = day_since_dose_1>=0 & day_since_dose_2 >=0 & day_since_dose_3 <0 & day_since_dose_4 <0;
VE_3_indicator = day_since_dose_1>=0 & day_since_dose_2 >=0 & day_since_dose_3 >=0 & day_since_dose_4 <0;
VE_4_indicator = day_since_dose_1>=0 & day_since_dose_2 >=0 & day_since_dose_3 >=0 & day_since_dose_4 >=0;

% Calculate weekly cumulative p_infection and immunity given VE estimates
if (no_dose==0)
    VE_mat = zeros(size(tv,1),1);
else
    VE_mat = VE_BNT1*exp(-waning_BNT*day_since_dose_1).*VE_1_indicator.*(type_dose_1==1)...
        +VE_C1*exp(-waning_C*day_since_dose_1).*VE_1_indicator.*(type_dose_1==2)...
        +VE_BNT2*exp(-waning_BNT*day_since_dose_2).*VE_2_indicator.*(type_dose_2==1)...
        +VE_C2*exp(-waning_C*day_since_dose_2).*VE_2_indicator.*(type_dose_2==2)...
        +VE_BNT3*exp(-waning_BNT*day_since_dose_3).*VE_3_indicator.*(type_dose_3==1)...
        +VE_C3*exp(-waning_C*day_since_dose_3).*VE_3_indicator.*(type_dose_3==2)...
        +VE_BNT4*exp(-waning_BNT*day_since_dose_4).*VE_4_indicator.*(type_dose_4==1)...
        +VE_C4*exp(-waning_C*day_since_dose_4).*VE_4_indicator.*(type_dose_4==2);
end

p_infection = zeros(1,numel(proxy));
immunity = zeros(1,numel(proxy));

proxy_vax = proxy.*(1-VE_mat);

cum_proxy_vax = cumsum(proxy_vax);
p_infection = 1-exp(-(age_scaling_factor.*total_FOI).*cum_proxy_vax);

immunity = p_infection + (1-p_infection).*VE_mat;

p_infection_week = p_infection(tv_week);
immunity_week = immunity(tv_week);

end