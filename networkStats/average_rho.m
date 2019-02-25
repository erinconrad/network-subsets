function rho_mean = average_rho(rho,dim)

%{
This function averages Spearman rank coefficients by way of doing Fisher's
transformation to get z values, then averaging the z values, and then
backtransforming to get an average rho
%}

z = atanh(rho);
z_mean = nanmean(z,dim);
rho_mean = tanh(z_mean);

end