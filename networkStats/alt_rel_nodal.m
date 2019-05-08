function reliability = alt_rel_nodal(perm_metric)

variance_error = nanmean(squeeze(nanstd(perm_metric,0,3).^2),1);

variance_true = nanmean((nanstd(perm_metric,0,1).^2),3);
% Take the variance across the electrodes for each random resampling and
% average that variance

reliability = variance_true./(variance_true + variance_error);


end 