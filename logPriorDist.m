function funval = logPriorDist(x,data)

VE_BNT2 = prod(x(2:4));
VE_C2 = prod(x(5:7));
age_factors = x(16:19);

half_life = log(2)./x(8:9);
funval = sum(-10*(max(0, half_life-365)+max(0, 7-half_life)));
funval = funval-100*max(0, VE_C2-VE_BNT2);
funval = funval-100*max(0, age_factors(1)-2)-100*max(0,age_factors(2)-2)-100*max(0,age_factors(3)-2)-100*max(0,age_factors(4)-2);



