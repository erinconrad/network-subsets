function clinical_corr(pt,stats)

outcome_all = [];
sync_std_80 = [];
sync_mean_80 = [];

%% Look at contiguous and -5 seconds
contig_text = 'random';
sec_text = 'sec_0';

for whichPt = 1:length(stats)
    
    % synchronizability
    sync_std_80_temp = stats(whichPt).sync.(contig_text).(sec_text).std(4);
    sync_std_80 = [sync_std_80;sync_std_80_temp]; 
    sync_mean_80 = [sync_mean_80;stats(whichPt).sync.(contig_text).(sec_text).mean(4)];
    
    % outcome
    outcome_all = [outcome_all;getOutcome(pt,whichPt)];
    
    
end

% First, make std relative to std among patients
std_patients_sync = std(sync_mean_80);
sync_std_rel = sync_std_80/std_patients_sync;

% correlate sync_std_80 with outcome
figure
scatter(sync_std_rel,outcome_all)
[rho,p] = corr(sync_std_rel,outcome_all,'Type','Spearman')

end