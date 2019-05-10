function reliability = alt_rel_nodal(perm_metric)
%{
THIS IS THE METHOD WE NOW USE TO CALCULATE NODAL RELIABILITY
%}

% This is the variance across random permutations, averaged over all
% electrodes
variance_error = nanmean(squeeze(nanstd(perm_metric,0,3).^2),1);

% This is the variance across electrodes for a given permutation of the
% graph at that removal percentage, then averaged over all permutations for
% that removal percentage
variance_true = nanmean((nanstd(perm_metric,0,1).^2),3);


reliability = variance_true./(variance_true + variance_error);


end 