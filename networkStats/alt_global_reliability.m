function reliability = alt_global_reliability(perm_metric)

%{
This is the code we now use for global reliability
%}

% This is the variance across patients for a given permutation within a
% removal percentage, then averaged over all permutations for that removal
% percentage. Then we repeat it so that every patient gets the same true
% variance (it is the same for everyone because since these are global
% metrics the true variance has to be defined across patients).
variance_true = repmat((nanmean(nanstd(perm_metric,0,1).^2,4)),...
    [size(perm_metric,1),1,1]);

% This is the variance of the global metric across permutations
variance_error = ((nanstd(perm_metric,0,4).^2));

reliability = variance_true./(variance_true+variance_error);

end


