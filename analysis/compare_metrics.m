function compare_metrics(stats)


[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'basic_metrics/'];

doPlots = 1;

%% Initialize arrays
nodal_metrics = {'cc','cc_reg','ns','bc','ec','clust'};
global_metrics = {'sync','eff','trans'};
ef = [20 40 60 80 100];

np = length(stats);

names ={};
name_nums = {};
ag_nodal = nan(np,length(nodal_metrics),length(ef));
var_nodal = nan(np,length(nodal_metrics),length(ef));
ag_global = nan(np,length(global_metrics),length(ef));
std_global = nan(np,length(global_metrics),length(ef));
true_global = nan(np,length(global_metrics));

for i = 1:length(stats)
    
    if isempty(stats(i).name) == 1, continue; end
    %% Extract just numbers from name (for plotting)
    names = [names;stats(i).name];
    [num_idx_s] = regexp(stats(i).name,'\d+');
    name_nums = [name_nums;stats(i).name(num_idx_s:end)];
    
    %% Look at contiguous and -5 seconds
    contig_text = 'random';
    sec_text = 'sec_neg5';
    freq = 'high_gamma';
    
    base = stats(i).(freq).(contig_text).(sec_text);
    
    %% Get agreement and variability for metrics by resection size
    for j = 1:length(nodal_metrics)
        var_nodal(i,j,:) = base.(nodal_metrics{j}).rel_std;
        ag_nodal(i,j,:) = base.(nodal_metrics{j}).rho_mean';
    end
    
    %% Get agreement and variability for global metrics by resection size
    for j = 1:length(global_metrics)
        std_global(i,j,:) = base.(global_metrics{j}).std';
        ag_global(i,j,:) = mean(base.(global_metrics{j}).rel_diff_norm,2);
        true_global(i,j,:) = base.(global_metrics{j}).true;
    end
   
end

%% Get var_global (relative std by dividing by std across patients)
var_global = std_global./nanstd(true_global,0,1);

%% Average over patients
avg_ag_nodal = squeeze(average_rho(ag_nodal,1)); % Fisher transform for rho
avg_var_nodal = squeeze(nanmean(var_nodal,1));
avg_ag_global = squeeze(nanmean(ag_global,1));
avg_var_global = squeeze(nanmean(var_global,1));

%% Individual patient, 80% retained
var_nodal_80 = var_nodal(:,:,4);
var_global_80 = var_global(:,:,4);

%% Calculate statistics for variability

% Compare variability when 80% retained for global metrics
[p,~,stats1] = signrank(var_global_80(:,1),var_global_80(:,2));
fprintf('Sign rank test for global metrics: p = %1.2e, sign-rank = %d\n',...
    p, stats1.signedrank);

% Compare variability when 80% retained for nodal metrics
[p,tbl] = friedman(var_nodal_80(~isnan(var_nodal_80(:,1)),[1,3:end]),1,'off');
fprintf('Friedman test for global metrics: p = %1.1e, chi-squared = %1.1f\n',...
    p, tbl{2,5});

[p,~,stats1] = signrank(var_nodal_80(:,1),var_nodal_80(:,2));
fprintf('Sign rank test for cc vs ns: p = %1.2e, sign-rank = %d\n',...
    p, stats1.signedrank);

[p,~,stats1] = signrank(var_nodal_80(:,1),var_nodal_80(:,3));
fprintf('Sign rank test for cc vs bc: p = %1.2e, sign-rank = %d\n',...
    p, stats1.signedrank);

[p,~,stats1] = signrank(var_nodal_80(:,2),var_nodal_80(:,3));
fprintf('Sign rank test for ns vs bc: p = %1.2e, sign-rank = %d\n',...
    p, stats1.signedrank);

%% Plot averages across patients
if doPlots == 1

    % Nodal metrics
    figure
    set(gcf,'Position',[174 207 1136 582])
    [ha, pos] = tight_subplot(2, 2, [0.07 0.1], [0.1 0.06],[0.11 0.02]);
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
    ylabel({'Spearman rank correlation',...
        'with original'})
    title('Average agreement by subsample size');
    set(gca,'Fontsize',20);

    % Variability
    axes(ha(2))
    for j = 1:size(avg_ag_nodal,1)
         scatter(ef,avg_var_nodal(j,:),200,'filled');
         hold on
    end
    legend('Control centrality','Regional control centrality',...
        'Node strength','Betweenness centrality',...
        'Eigenvector centrality','Clustering coefficient',...
        'location','northeast');
    %xlabel('Percent nodes retained');
    ylabel({'Relative standard deviation'})
    title('Variability by subsample size');
    set(gca,'Fontsize',20);


    % Global metrics
    % Average agreement
    axes(ha(3))
    for j = 1:size(avg_ag_global,1)
         scatter(ef,avg_ag_global(j,:),200,'filled');
         hold on
    end
    legend('Synchronizability','Global efficiency','Transitivity',...
        'location','northeast');
    xlabel('Percent nodes retained');
    ylabel({'Relative difference',...
        'from original'})
    %title('Average agreement by subsample size');
    set(gca,'Fontsize',20);


    % Variability
    axes(ha(4))
    for j = 1:size(avg_ag_global,1)
         scatter(ef,avg_var_global(j,:),200,'filled');
         hold on
    end
    legend('Synchronizability','Global efficiency','Transitivity',...
        'location','northeast');
    xlabel('Percent nodes retained');
    ylabel({'Relative standard deviation'})
    %title('Variability by subsample size');
    set(gca,'Fontsize',20);
    pause
    print(gcf,[outFolder,'avg_metrics_',freq,contig_text,sec_text],'-depsc');
    close(gcf)

    %% Plot 80% for all patients

    figure
    set(gcf,'Position',[31 173 1400 632]);
    [ha, pos] = tight_subplot(2, 1, [0.1 0.15], [0.1 0.07],[0.06 0.02]);
    cols = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
        0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330];

    % Nodal metrics
    axes(ha(1))
    count = 0;
    for i = 1:size(var_nodal_80,1)
        if isempty(stats(i).name) == 1, continue; end
        count = count + 1;
        for j = 1:size(var_nodal_80,2)
            scatter(count,var_nodal_80(i,j),200,cols(j,:),'filled');
            hold on
        end
    end
    legend('Control centrality','Regional control centrality',...
        'Node strength','Betweenness centrality',...
        'Eigenvector centrality','Clustering coefficient');
    xticks(1:length(names))
    %xlabel('Which patient')
    ylabel('Relative standard deviation')
    set(gca,'fontsize',20)
    title('Variability of metric when 20% of network removed');
    %xticklabels(name_nums)


    % Global metrics
    axes(ha(2))
    count = 0;
    for i = 1:size(var_global_80,1)
        if isempty(stats(i).name) == 1, continue; end
        count = count + 1;
        for j = 1:size(var_global_80,2)
            scatter(count,var_global_80(i,j),200,cols(j,:),'filled');
            hold on
        end
    end
    legend('Synchronizability','Global efficiency','Transitivity');
    xticks(1:length(names))
    xlabel('Which patient')
    ylabel('Relative standard deviation')
    set(gca,'fontsize',20)
    pause
    print(gcf,[outFolder,'all_pat_80_',freq,contig_text,sec_text],'-depsc');
    close(gcf)

end

%xticklabels(name_nums)


end