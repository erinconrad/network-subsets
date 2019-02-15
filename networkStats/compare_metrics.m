function compare_metrics(stats)


%% Initialize arrays
names = {};
sync_std_80 = [];
eff_std_80 = [];
cc_rho_80 = [];
ns_rho_80 = [];
bc_rho_80 = [];

for i = 1:length(stats)
    
    names = [names;stats(i).name];
    
    %% Look at contiguous and -5 seconds
    contig_text = 'contig';
    sec_text = 'sec_neg5';
    
    %% For global metrics, get std at 80%
    
    % synchronizability
    sync_std_80_temp = stats(i).sync.(contig_text).(sec_text).std(4);
    sync_std_80 = [sync_std_80;sync_std_80_temp];
    
    % global efficiency
    eff_std_80_temp = stats(i).eff.(contig_text).(sec_text).std(4);
    eff_std_80 = [eff_std_80;eff_std_80_temp];
    
    %% For node metrics, get rho at 80%
    
    % control centrality
    cc_rho_80_temp = stats(i).cc.(contig_text).(sec_text).rho.mean(4);
    cc_rho_80 = [cc_rho_80;cc_rho_80_temp];
    
    % node strength
    ns_rho_80_temp = stats(i).ns.(contig_text).(sec_text).rho.mean(4);
    ns_rho_80 = [ns_rho_80;ns_rho_80_temp];
    
    % betweenness centrality
    bc_rho_80_temp = stats(i).bc.(contig_text).(sec_text).rho.mean(4);
    bc_rho_80 = [bc_rho_80;bc_rho_80_temp];
    
end

%% Plot these

% Global metrics
figure
subplot(1,2,1)
plot(sync_std_80)
xticklabels(names)
subplot(1,2,2)
plot(eff_std_80)
xticklabels(names)

% Node metrics
figure
subplot(1,3,1)
plot(cc_rho_80)
xticklabels(names)
subplot(1,3,2)
plot(ns_rho_80)
xticklabels(names)
subplot(1,3,3)
plot(bc_rho_80)
xticklabels(names)



end