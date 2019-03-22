function reliability_nodal(perm_metric,true_metric)

%{
Total variance = variance of true scores + variance of error
Reliability =  variance of true scores/total variance
%}

variance_error = nanmean(squeeze(nanstd(perm_metric,0,3).^2),1);

variance_true = nanstd(true_metric,0,1).^2;

reliability = variance_true/(variance_true + variance_total);


end