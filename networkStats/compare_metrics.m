function compare_metrics(stats)


[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'basic_metrics/'];

%% Initialize arrays
names = {};
name_nums = {};
sync_std_80 = [];
eff_std_80 = [];
sync_mean_80 = [];
eff_mean_80 = [];
cc_rho_80 = [];
ns_rho_80 = [];
bc_rho_80 = [];
resect_wrong = [];
sync_std_all = [];
eff_std_all = [];
cc_rho_all = [];
ns_rho_all = [];
bc_rho_all = [];
sync_mean_all = [];
eff_mean_all = [];

for i = 1:length(stats)
    
    %% Extract just numbers from name (for plotting)
    names = [names;stats(i).name];
    [num_idx_s] = regexp(stats(i).name,'\d+');
    name_nums = [name_nums;stats(i).name(num_idx_s:end)];
    
    %% Look at contiguous and -5 seconds
    contig_text = 'contiguous';
    sec_text = 'sec_neg5';
    
    
    %% For global metrics, get std at 80%
    
    % synchronizability
    sync_std_80_temp = stats(i).sync.(contig_text).(sec_text).std(4);
    sync_std_80 = [sync_std_80;sync_std_80_temp];
    sync_std_all = [sync_std_all,stats(i).sync.(contig_text).(sec_text).std];
    
    % global efficiency
    eff_std_80_temp = stats(i).eff.(contig_text).(sec_text).std(4);
    eff_std_80 = [eff_std_80;eff_std_80_temp];
    eff_std_all = [eff_std_all,stats(i).eff.(contig_text).(sec_text).std];
    
    sync_mean_80 = [sync_mean_80;stats(i).sync.(contig_text).(sec_text).mean(4)];
    eff_mean_80 = [eff_mean_80;stats(i).eff.(contig_text).(sec_text).mean(4)];
    sync_mean_all = [sync_mean_all,stats(i).sync.(contig_text).(sec_text).mean];
    eff_mean_all = [eff_mean_all,stats(i).eff.(contig_text).(sec_text).mean];
    
    
    %% For node metrics, get rho at 80%
    
    % control centrality
    cc_rho_80_temp = stats(i).cc.(contig_text).(sec_text).rho.mean(4);
    cc_rho_80 = [cc_rho_80;cc_rho_80_temp];
    cc_rho_all = [cc_rho_all,stats(i).cc.(contig_text).(sec_text).rho.mean];
    
    % node strength
    ns_rho_80_temp = stats(i).ns.(contig_text).(sec_text).rho.mean(4);
    ns_rho_80 = [ns_rho_80;ns_rho_80_temp];
    ns_rho_all = [ns_rho_all,stats(i).ns.(contig_text).(sec_text).rho.mean];
    
    % betweenness centrality
    bc_rho_80_temp = stats(i).bc.(contig_text).(sec_text).rho.mean(4);
    bc_rho_80 = [bc_rho_80;bc_rho_80_temp];
    bc_rho_all = [bc_rho_all,stats(i).bc.(contig_text).(sec_text).rho.mean];
    
    % How often do we resect wrong brain if keep 80% electrodes
    resect_wrong = [resect_wrong;stats(i).cc.(contig_text).(sec_text).resect_wrong(4)];
    
   
    
end

%% Do statistics to see if one is better than others

%% Global metrics: I can use a Wilcoxon signed-rank test since just 2
% I don't know if it's fair to compare standard deviations of these 2
% different measures

% First, make std relative to std among patients
std_patients_sync = std(sync_mean_80);
std_patients_eff = std(eff_mean_80);
sync_std_rel = sync_std_80/std_patients_sync;
eff_std_rel = eff_std_80/std_patients_eff;
sync_std_all_rel = sync_std_all./std(sync_mean_all,0,2);
eff_std_all_rel = eff_std_all./std(eff_mean_all,0,2);

[p,~,stats0] = signrank(sync_std_80,eff_std_80);
fprintf(['Average std of metric when 80%% of network retained is\n'...
    '%1.2f for synchronizability and %1.2f for global efficiency.\n'...
    'Wilcoxon signed rank: p = %1.2e and sign rank = %1.1f\n\n'],...
    mean(sync_std_80), mean(eff_std_80),p,stats0.signedrank);

[p,~,stats0] = signrank(sync_std_rel,eff_std_rel);
fprintf(['Average relative std of metric when 80%% of network retained is\n'...
    '%1.2f for synchronizability and %1.2f for global efficiency.\n'...
    'Wilcoxon signed rank: p = %1.2e and sign rank = %1.1f\n\n'],...
    mean(sync_std_rel), mean(eff_std_rel),p,stats0.signedrank);

%{
[p,~,stats0] = signrank(sync_mean_80,eff_mean_80);
fprintf(['Average mean of metric when 80%% of network retained is\n'...
    '%1.2f for synchronizability and %1.2f for global efficiency.\n'...
    'Wilcoxon signed rank: p = %1.2e and sign rank = %1.1f\n\n'],...
    mean(sync_mean_80), mean(eff_mean_80),p,stats0.signedrank);
    %}

%% Node metrics: 3, so need Friedman test
[p,tbl] = friedman([cc_rho_80,ns_rho_80,bc_rho_80],1,'off');
[p_cc_ns,~,stats0] = signrank(ns_rho_80,cc_rho_80);
[p_cc_bc,~,stats1] = signrank(bc_rho_80,cc_rho_80);
[p_ns_bc,~,stats2] = signrank(ns_rho_80,bc_rho_80);
fprintf(['Average rho when 80%% of network retained is\n'...
    '%1.2f for control centrality, %1.2f for node strength\n,'...
    'and %1.2f for betweenness centrality.\n'...
    'Friedman test: p = %1.2e and Q = %1.2f\n\n'...
    'CC-NS signed-rank test: p = %1.2e and sign rank = %1.1f\n'...
    'CC-BC signed-rank test: p = %1.2e and sign rank = %1.1f\n'...
    'NS-BC signed-rank test: p = %1.2e and sign rank = %1.1f\n\n'],...
    mean(cc_rho_80), mean(ns_rho_80), mean(bc_rho_80),...
    p,tbl{2,5},p_cc_ns,stats0.signedrank,p_cc_bc,stats1.signedrank,p_ns_bc,stats2.signedrank);

%% All plots, average across patients

% Global
figure
set(gcf,'Position',[1 432 800 366]);
plot(1:size(sync_std_all_rel,1),mean(sync_std_all_rel,2),'LineWidth',2)
hold on
plot(1:size(eff_std_all_rel,1),mean(eff_std_all_rel,2),'LineWidth',2)
legend({'Synchronizability','Global efficiency'},'location','southeast')
xticks(1:size(sync_std_all_rel,1))
ylabel('Standard deviation')
xlabel('Percent of network retained');
xticklabels({'20%','40%','60%','80%','100%'})
title(sprintf(['Relative standard deviation of metric, averaged across all patients\n'...
    '%s, %s'],contig_text,sec_text),'Interpreter','none');
set(gca,'fontsize',15)
print([outFolder,sprintf('global_allperc_%s_%s',contig_text,sec_text)],'-depsc')

% Nodal
figure
set(gcf,'Position',[1 432 800 366]);
plot(1:size(cc_rho_all,1),mean(cc_rho_all,2),'LineWidth',2)
hold on
plot(1:size(ns_rho_all,1),mean(ns_rho_all,2),'LineWidth',2)
plot(1:size(bc_rho_all,1),mean(bc_rho_all,2),'LineWidth',2)
legend({'Control centrality','Node strength','Betweenness centrality'},'location','southeast')
xticks(1:size(sync_std_all_rel,1))
xticklabels({'20%','40%','60%','80%','100%'})
xlabel('Percent of network retained');
ylabel('Spearman rank coefficient')
title(sprintf(['Spearman rank coefficient across electrodes between original metric\n'...
    'and metric, averaged across all patients\n%s, %s'],contig_text,sec_text),...
    'Interpreter','none');
set(gca,'fontsize',15)
print([outFolder,sprintf('nodal_allperc_%s_%s',contig_text,sec_text)],'-depsc')

%% 80% plots
% Global metrics
figure
set(gcf,'Position',[1 432 1440 366]);
scatter(1:length(names),sync_std_rel,100,'filled')
hold on
scatter(1:length(names),eff_std_rel,100,'filled')
legend({'Synchronizability','Global efficiency'},'location','southeast')
xticks(1:length(names))
xticklabels(name_nums)
title(sprintf(['Relative standard deviation of metric when 80%% of electrodes retained\n'...
    '%s, %s'],contig_text,sec_text),'Interpreter','none');
xlabel('Which patient')
ylabel('Standard deviation')
set(gca,'fontsize',15)
print([outFolder,sprintf('global_%s_%s',contig_text,sec_text)],'-depsc')


% Node metrics
figure
set(gcf,'Position',[1 432 1440 366]);
scatter(1:length(names),cc_rho_80,100,'filled')
hold on
scatter(1:length(names),ns_rho_80,100,'filled')
scatter(1:length(names),bc_rho_80,100,'filled')
legend({'Control centrality','Node strength','Betweenness centrality'},'location','southeast')
xticks(1:length(names))
xticklabels(name_nums)
title(sprintf(['Spearman rank coefficient across electrodes between original metric\n'...
    'and metric when 80%% of electrodes retained\n%s, %s'],contig_text,sec_text),...
    'Interpreter','none');
xlabel('Which patient')
ylabel('Spearman rank coefficient')
set(gca,'fontsize',15)
print([outFolder,sprintf('node_%s_%s',contig_text,sec_text)],'-depsc')

% How often do we resect wrong brain?
figure
set(gcf,'Position',[1 432 1440 366]);
scatter(1:length(names),resect_wrong,100,'filled')
xticks(1:length(names))
xticklabels(name_nums)
title(sprintf(['%% of time a desynchronizing node is labeled the most synchronizing\n'...
    'when 80%% of electrodes retained\n%s, %s'],contig_text,sec_text),...
    'Interpreter','none');
xlabel('Which patient')
ylabel('% of time')
set(gca,'fontsize',15)
print([outFolder,sprintf('resect_wrong_%s_%s',contig_text,sec_text)],'-depsc')





end