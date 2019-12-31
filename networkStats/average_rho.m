function [rho_mean,rho_var,z_mean,z_var] = average_rho(rho,dim)

%{
This function averages Spearman rank coefficients by way of doing Fisher's
transformation to get z values, then averaging the z values, and then
backtransforming to get an average rho
%}

z = atanh(rho);
z_mean = nanmean(z,dim);
z_var = (nanstd(z,0,dim)).^2;
rho_mean = tanh(z_mean);
rho_var = tanh(z_var);

end