function compare_metrics(pt,stats)

%{
This function takes the output of the resampling script network_stats and
compares network measure reliability
%}

%% Parameters
doPlots = 1; % plot things?
all_contig = {'random','contiguous'}; % look at random or contiguous removal
all_freq = {'high_gamma','beta'}; % which frequency coherence
all_sec = {'sec_neg10','sec_neg5','sec_0','sec_5','sec_10'}; % which times relative to EEC



%% Locations

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'compare_metrics/'];
addpath(genpath(scriptFolder))


%% Initialize arrays
nodal_metrics = {'cc','ns','bc','ec','clust'};
global_metrics = {'sync','eff','trans'};
ef = [20 40 60 80 100];

np = length(stats);

%% Initialize cross time, freq, etc. comparisons
order_global = nan(length(global_metrics),length(all_sec),...
    length(all_freq),length(all_contig));
order_nodal = nan(length(nodal_metrics),length(all_sec),...
    length(all_freq),length(all_contig));

all_global = nan(length(global_metrics),length(all_sec),...
    length(all_freq),length(all_contig));
all_nodal = nan(length(nodal_metrics),length(all_sec),...
    length(all_freq),length(all_contig));

% Loop through contig vs random removal, frequencies, and times
for contig_idx = 1%1:length(all_contig)
for freq_idx = 1%1:length(all_freq)
for sec_idx = 3%1:length(all_sec)
    
% Get appropriate contig vs random, frequency, time
contig_text = all_contig{contig_idx};
sec_text = all_sec{sec_idx};
freq = all_freq{freq_idx};
    
% Initialize arrays
names ={};
name_nums = {};
ag_nodal = nan(np,length(nodal_metrics),length(ef));
var_nodal = nan(np,length(nodal_metrics),length(ef));
ag_global = nan(np,length(global_metrics),length(ef));
std_global = nan(np,length(global_metrics),length(ef));
global_all = nan(np,length(global_metrics),length(ef),1000);
true_global = nan(np,length(global_metrics));
n_elecs = nan(np,1);

same_hub_nodal = nan(np,length(nodal_metrics),length(ef));

% Loop through patients
for i = 1:length(stats)
    
    if isempty(stats(i).name) == 1
        if doPlots == 0
            names = [names;nan]; 
        end
        continue; 
    end
    
    %% Get number of electrodes
    if strcmp(pt(i).name,stats(i).name) == 0, error('what\n'); end
    n_elecs(i) = length(pt(i).new_elecs.electrodes);
    
    %% Extract just numbers from name (for plotting)
    names = [names;stats(i).name];
    [num_idx_s] = regexp(stats(i).name,'\d+');
    name_nums = [name_nums;stats(i).name(num_idx_s:end)];

    if isfield(stats(i).(freq).(contig_text),sec_text) == 0, continue; end
    base = stats(i).(freq).(contig_text).(sec_text);
    
    %% Get agreement and reliability for metrics by resection size
    for j = 1:length(nodal_metrics)
        var_nodal(i,j,:) = base.(nodal_metrics{j}).rel_alt;
        ag_nodal(i,j,:) = base.(nodal_metrics{j}).rho_mean';
        
        %same_hub_nodal(i,j,:) = mean(base.(nodal_metrics{j}).same_hub,2);
    end
    
    %% Get agreement and variability for global metrics by resection size
    for j = 1:length(global_metrics)
        std_global(i,j,:) = base.(global_metrics{j}).std';
        ag_global(i,j,:) = mean(base.(global_metrics{j}).rel_diff_norm,2);
        true_global(i,j,:) = base.(global_metrics{j}).true;
        global_all(i,j,1:size(base.(global_metrics{j}).all,1),...
            1:size(base.(global_metrics{j}).all,2)) = base.(global_metrics{j}).all;
    end
   
end

%% Get reliability for global metrics
% nanstd(true_global,0,1) is the standard deviation of the global metric
% across patients
%var_global = global_reliability(std_global,nanstd(true_global,0,1));
var_global = alt_global_reliability(global_all);

%% Average over patients
avg_ag_nodal = squeeze(average_rho(ag_nodal,1)); % Fisher transform for rho
avg_var_nodal = squeeze(nanmean(var_nodal,1));
avg_ag_global = squeeze(nanmean(ag_global,1));
avg_var_global = squeeze(nanmean(var_global,1));
%avg_samehub_nodal = squeeze(nanmean(same_hub_nodal,1));

%% Individual patient, 80% retained
% These are the reliability metrics
var_nodal_80 = var_nodal(:,:,4);
var_global_80 = var_global(:,:,4);
%same_hub_nodal_80 = same_hub_nodal(:,:,4);

%% Calculate statistics for variability
if sum(sum(isnan(var_global_80))) == sum(sum((ones(size(var_global_80)))))
    continue;
end

% Descriptive stats for global metrics, 80% retained
for i = 1:length(global_metrics)
    fprintf(['The mean reliability when retaining 80%% of network for %s\n'...
        'is %1.3f\n\n'],global_metrics{i},nanmean(var_global_80(:,i)));
end

% Compare variability when 80% retained for global metrics
[p,tbl,stats1] = friedman(var_global_80(~isnan(var_global_80(:,1)),:),1,'off');
fprintf('Friedman test for global metrics: p = %1.1e, chi-squared = %1.1f, dof = %d\n',...
    p, tbl{2,5},tbl{2,3});

% perform a post-hoc Dunn's test
[c,~,~,~,t] = multcompare_erin(stats1,'CType','dunn-sidak','Display','off');
for i = 1:size(c,1)
    fprintf('Dunn''s test comparing %s and %s: t = %1.2f, p = %1.3f\n\n',...
        global_metrics{c(i,1)},global_metrics{c(i,2)},t(i),c(i,6));   
end

% compare this do doing a signed rank
[p,~,stats1] = signrank(var_global_80(:,1),var_global_80(:,2));
fprintf('Sign rank test for sync vs eff: p = %1.2e, sign-rank = %d\n',...
    p, stats1.signedrank);

[p,~,stats1] = signrank(var_global_80(:,1),var_global_80(:,3));
fprintf('Sign rank test for sync vs transitivity: p = %1.2e, sign-rank = %d\n',...
    p, stats1.signedrank);

[p,~,stats1] = signrank(var_global_80(:,2),var_global_80(:,3));
fprintf('Sign rank test for eff vs transitivity: p = %1.2e, sign-rank = %d\n',...
    p, stats1.signedrank);

% Descriptive stats for nodal metrics, 80% retained
for i = 1:length(nodal_metrics)
    fprintf(['The mean reliability when retaining 80%% of network for %s\n'...
        'is %1.2f\n\n'],nodal_metrics{i},nanmean(var_nodal_80(:,i)));
end

% Compare variability when 80% retained for nodal metrics
[p,tbl,stats1] = friedman(var_nodal_80(~isnan(var_nodal_80(:,1)),:),1,'off');
fprintf('Friedman test for nodal metrics: p = %1.1e, chi-squared = %1.1f, dof = %d\n',...
    p, tbl{2,5},tbl{3,3});

% perform a post-hoc Dunn's test
[c,~,~,~,t] = multcompare_erin(stats1,'CType','dunn-sidak','Display','off');
for i = 1:size(c,1)
    fprintf('Dunn''s test comparing %s and %s: t = %1.2f, p = %1.3f\n\n',...
        nodal_metrics{c(i,1)},nodal_metrics{c(i,2)},t(i),c(i,6));   
end

% Descriptive stats for how often hub stays the same
%{
for i = 1:length(nodal_metrics)
    fprintf(['The hub stability when retaining 80%% of the network for %s\n'...
        'is %1.2f\n\n'],nodal_metrics{i},nanmean(same_hub_nodal_80(:,i)));    
end
%}

% Compare hub stability when 80% retained for nodal metrics
%{
[p,tbl,stats1] = friedman(same_hub_nodal_80(~isnan(same_hub_nodal_80(:,1)),:),1,'off');
fprintf('Friedman test for nodal hub stability: p = %1.1e, chi-squared = %1.1f, dof = %d\n',...
    p, tbl{2,5},tbl{3,3});
%}

% perform a post-hoc Dunn's test
[c,~,~,~,t] = multcompare_erin(stats1,'CType','dunn-sidak','Display','off');
for i = 1:size(c,1)
    fprintf('Dunn''s test comparing hub stability for %s and %s: t = %1.2f, p = %1.3f\n\n',...
        nodal_metrics{c(i,1)},nodal_metrics{c(i,2)},t(i),c(i,6));   
end


%% Order the nodal and global metrics
% Calculate mean reliability
mean_rel_global = mean(var_global_80(~isnan(var_global_80(:,1)),:),1);
mean_rel_nodal = mean(var_nodal_80(~isnan(var_nodal_80(:,1)),:),1);

[~,order_global(:,sec_idx,freq_idx,contig_idx)] = sort(mean_rel_global,'descend'); 
[~,order_nodal(:,sec_idx,freq_idx,contig_idx)] = sort(mean_rel_nodal,'descend'); 

all_global(:,sec_idx,freq_idx,contig_idx) = mean_rel_global;
all_nodal(:,sec_idx,freq_idx,contig_idx) = mean_rel_nodal;

mean_rel_global_all_removal_perc = squeeze(nanmean(var_global,1));
mean_rel_nodal_all_removal_perc = squeeze(nanmean(var_nodal,1));

%% Correlate reliability with number of electrodes
% global reliability
rho_global = zeros(3,1);
for i = 1:size(var_global_80,2)
    [rho,pval] = corr(var_global_80(~isnan(var_global_80(:,1)),i),n_elecs(~isnan(var_global_80(:,1))),'Type','Spearman'); 
    
    fprintf(['Spearman rank correlation between %s reliability and'...
        'electrode number:\n rho = %1.2f, df = %d, p = %1.3f\n.'],...
        global_metrics{i},rho,length(n_elecs(~isnan(var_global_80(:,1))))-2,...
        pval);
    
    rho_global(i) = rho;
    
    if 1==0
        scatter(var_global_80(:,i),n_elecs)
        text(0.9,50,sprintf('%1.2f',rho),'fontsize',20)
        pause
        close(gcf)
    end
end
rho_global_avg = average_rho(rho_global,1);
fprintf(['The correlation coefficient between global reliability and\n'...
    'electrode number is %1.2f\n'],rho_global_avg);

% nodal reliability
rho_nodal = zeros(5,1);
for i = 1:size(var_nodal_80,2)
    [rho,pval] = corr(var_nodal_80(~isnan(var_nodal_80(:,1)),i),n_elecs(~isnan(var_nodal_80(:,1))),'Type','Spearman'); 
    rho_nodal(i) = rho;
    fprintf(['Spearman rank correlation between %s reliability and'...
        'electrode number:\n rho = %1.2f, df = %d, p = %1.3f\n.'],...
        nodal_metrics{i},rho,length(n_elecs(~isnan(var_nodal_80(:,1))))-2,...
        pval);
    if 1==0
        scatter(var_nodal_80(:,i),n_elecs)
        text(0.9,50,sprintf('%1.2f',rho),'fontsize',20)
        pause
        close(gcf)
    end
end
rho_nodal_avg = average_rho(rho_nodal,1);
fprintf(['The correlation coefficient between nodal reliability and\n'...
    'electrode number is %1.2f\n'],rho_nodal_avg);


%% Reliability for multiple removal percentages
mean_rel_nodal_all_removal_perc
mean_rel_global_all_removal_perc


end
end
end

%% Determine consistency in order
% Keep random and high gamma, vary time
order_global(:,:,1,1)
order_nodal(:,:,1,1)

% Keep time and random, vary freq
squeeze(order_global(:,1,:,1))
squeeze(order_nodal(:,1,:,1))

% Keep time and freq, vary random vs contig
squeeze(order_global(:,1,1,:))
squeeze(order_nodal(:,1,1,:))

%% Tables for multiple times
fprintf('Multiple times:\n');
all_global(:,:,1,1)
all_nodal(:,:,1,1)

%% Tables for beta band
fprintf('EEC, beta band:\n');
all_nodal(:,3,2,1)
all_global(:,3,2,1)


%% Tables for contig removal
fprintf('EEC, contig removal:\n')
all_global(:,3,1,2)
all_nodal(:,3,1,2)

%% Tables for EEC, high gamma, random (sz 2)
fprintf('EEC, high gamma, random:\n');
all_global(:,3,1,1)
all_nodal(:,3,1,1)

%% Plot averages across patients
if doPlots == 1

    %{
    % Nodal metrics
    figure
    set(gcf,'Position',[174 207 1300 02])
    pos_f = get(gcf,'Position');
    [ha, pos] = tight_subplot(2, 2, [0.07 0.07], [0.1 0.06],[0.11 0.02]);
    % Average agreement
    axes(ha(1))
    for j = 1:size(avg_ag_nodal,1)
         scatter(ef,avg_ag_nodal(j,:),200,'filled');
         hold on
    end
    
    legend('Control centrality','Regional control centrality',...
        'Node strength','Betweenness centrality',...
        'Eigenvector centrality','Clustering coefficient',...
        'location','southeast');
    
    %xlabel('Percent nodes retained');
    ylabel({'NODAL METRICS','','Spearman rank correlation',...
        'with original'})
    title('Average agreement by subsample size');
    set(gca,'Fontsize',20);
    


    % Variability
    axes(ha(2))
    nd = zeros(size(avg_ag_nodal,1),1);
    for j = 1:size(avg_ag_nodal,1)
       nd(j) = scatter(ef,avg_var_nodal(j,:),200,'filled');
         hold on
    end
    %{
    legend('Control centrality','Regional control centrality',...
        'Node strength','Betweenness centrality',...
        'Eigenvector centrality','Clustering coefficient',...
        'location','northeastoutside');
    %}
    %xlabel('Percent nodes retained');
    ylabel({'Reliability'})
    title('Reliability by subsample size');
    set(gca,'Fontsize',20);
    
   % axes(ha(3))
   l1 = legend(nd,{'Control centrality',...
        'Node strength','Betweenness centrality',...
        'Eigenvector centrality','Clustering coefficient'},'Location',...
        'southeast');
   
    legend boxoff
   % set(gca,'Fontsize',20);

    % Global metrics
    % Average agreement
    axes(ha(3))
    for j = 1:size(avg_ag_global,1)
         scatter(ef,avg_ag_global(j,:),200,'filled');
         hold on
    end
    %{
    legend('Synchronizability','Global efficiency','Transitivity',...
        'location','northeast');
    %}
    xlabel('Percent nodes retained');
    ylabel({'GLOBAL METRICS','','Relative difference',...
        'from original'})
    %title('Average agreement by subsample size');
    set(gca,'Fontsize',20);


    % Variability
    gl = zeros(size(avg_ag_global,1),1);
    axes(ha(4))
    for j = 1:size(avg_ag_global,1)
         gl(j) = scatter(ef,avg_var_global(j,:),200,'filled');
         hold on
    end
    %{
    legend('Synchronizability','Global efficiency','Transitivity',...
        'location','northeastoutside');
    %}
    xlabel('Percent nodes retained');
    ylabel({'Reliability'})
    %title('Variability by subsample size');
    set(gca,'Fontsize',20);
    
   % axes(ha(6))
   l2= legend(gl,{'Synchronizability','Global efficiency','Transitivity'},'Location',...
        'southeast');
    legend boxoff
    %set(gca,'Fontsize',20);
    for i = 1:length(l1.EntryContainer.NodeChildren)
        l1.EntryContainer.NodeChildren(i).Icon.Transform.Children.Children.Size = 14;
    end
    
    % Necessary voodoo
    pause(3)
    
    for i = 1:length(l2.EntryContainer.NodeChildren)
        l2.EntryContainer.NodeChildren(i).Icon.Transform.Children.Children.Size = 14;
    end
    pause
    print(gcf,[outFolder,'avg_metrics_',freq,contig_text,sec_text],'-depsc');
    close(gcf)
    %}
    
    %% Just reliability
    
    cols = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
        0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;...
        0.6350 0.0780 0.1840;0.75 0.5 0.5];
    
    figure
    set(gcf,'Position',[174 207 1300 660])
    [ha, pos] = tight_subplot(2, 1, [0.11 0.04], [0.11 0.06],[0.06 0.02]);
    
    
    % Nodal reliability
    axes(ha(1))
    nd = zeros(size(avg_ag_nodal,1),1);
    for j = 1:size(avg_ag_nodal,1)
        std_rel = squeeze(nanstd(var_nodal,0,1));
       errorbar(100-ef+j*2-6,avg_var_nodal(j,:),std_rel(j,:),...
           'o','MarkerSize',15,'MarkerEdgeColor',cols(j,:),...
           'MarkerFaceColor',cols(j,:),'Color',cols(j,:),'linewidth',2);
         hold on
        nd(j) = plot(100-ef+j*2-6,avg_var_nodal(j,:),'o','MarkerSize',15,...
            'MarkerEdgeColor',cols(j,:),...
           'MarkerFaceColor',cols(j,:));
    end
    %{
    legend('Control centrality','Regional control centrality',...
        'Node strength','Betweenness centrality',...
        'Eigenvector centrality','Clustering coefficient',...
        'location','northeastoutside');
    %}
   % xlabel('Percent nodes removed');
    ylabel({'Reliability'})
    title('Nodal reliability');
    xlim([-6 86])
    xticks(sort(100-ef))
   % title('Reliability by subsample size','Position',[0.1 0.1 0.1 0.1]);
    set(gca,'Fontsize',25);
    
   % axes(ha(3))
   legend(nd,{'Control centrality',...
        'Node strength','Betweenness centrality',...
        'Eigenvector centrality','Clustering coefficient'},'Location',...
        'southwest','fontsize',25);
    legend boxoff
    
    gl = zeros(size(avg_ag_global,1),1);
    axes(ha(2))
    for j = 1:size(avg_ag_global,1)
         std_rel = squeeze(nanstd(var_global,1));
         errorbar(100-ef+j*2-4,avg_var_global(j,:),std_rel(j,:),...
           'o','MarkerSize',15,'MarkerEdgeColor',cols(j+5,:),...
           'MarkerFaceColor',cols(j+5,:),'Color',cols(j+5,:),'linewidth',2);
         hold on
        gl(j) = plot(100-ef+j*2-4,avg_var_global(j,:),'o','MarkerSize',15,...
            'MarkerEdgeColor',cols(j+5,:),...
           'MarkerFaceColor',cols(j+5,:));
    end
    %{
    legend('Synchronizability','Global efficiency','Transitivity',...
        'location','northeastoutside');
    %}
    xlabel('Percent nodes removed');
    ylabel({'Reliability'})
    title('Global reliability');
    xlim([-6 86])
    xticks(sort(100-ef))
    %title('Variability by subsample size');
    set(gca,'Fontsize',25);
    
   % axes(ha(6))
   legend(gl,{'Synchronizability','Global efficiency','Transitivity'},'Location',...
        'southwest','fontsize',25);
    legend boxoff
    
    annotation('textbox',[0 0.89 0.1 0.1],'String',...
        'A','FontSize',35,'linestyle','none');
    annotation('textbox',[0 0.40 0.1 0.1],'String',...
        'B','FontSize',35,'linestyle','none');
    
    
    %pause
    %print(gcf,[outFolder,'new_all_fig2_',freq,contig_text,sec_text],'-depsc');
    %close(gcf)

end

%xticklabels(name_nums)


end