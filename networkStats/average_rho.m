function all_avg = average_rho(rho,dim)

%{
This function averages Spearman rank coefficients by way of doing Fisher's
transformation to get z values, then averaging the z values, and then
backtransforming to get an average rho
%}

n = size(rho,dim); 
z = atanh(rho);

% Do a fix for case when rho = 1 (which would result in z == Inf without
% correction). The probability of this occurring between two random equal
% sized vectors due to chance is 1/n!. And so I will set z to be the
% z-value corresponding to a p-value of 1/n!. 
z_inf = -norminv(1/factorial(n));

z(z==inf) = z_inf;

z_mean = nanmean(z,dim);
z_var = (nanstd(z,0,dim)).^2;
rho_mean = tanh(z_mean);
rho_var = tanh(z_var);

all_avg.fisher_rho = rho_mean;
all_avg.fisher_var = rho_var;
all_avg.simple_rho = nanmean(rho,dim);
all_avg.simple_var = (nanstd(rho,0,dim)).^2;

end