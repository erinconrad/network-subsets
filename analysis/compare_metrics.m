function compare_metrics(stats)

%% Parameters
contig_text = 'random';
sec_text = 'sec_neg5';
freq = 'high_gamma';
doPlots = 1;

%% Locations

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'compare_metrics/'];


%% Initialize arrays
nodal_metrics = {'cc','ns','bc','ec','clust'};
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

% Loop through patients
for i = 1:length(stats)
    
    if isempty(stats(i).name) == 1, continue; end
    %% Extract just numbers from name (for plotting)
    names = [names;stats(i).name];
    [num_idx_s] = regexp(stats(i).name,'\d+');
    name_nums = [name_nums;stats(i).name(num_idx_s:end)];

    base = stats(i).(freq).(contig_text).(sec_text);
    
    %% Get agreement and reliability for metrics by resection size
    for j = 1:length(nodal_metrics)
        var_nodal(i,j,:) = base.(nodal_metrics{j}).rel;
        ag_nodal(i,j,:) = base.(nodal_metrics{j}).rho_mean';
    end
    
    %% Get agreement and variability for global metrics by resection size
    for j = 1:length(global_metrics)
        std_global(i,j,:) = base.(global_metrics{j}).std';
        ag_global(i,j,:) = mean(base.(global_metrics{j}).rel_diff_norm,2);
        true_global(i,j,:) = base.(global_metrics{j}).true;
    end
   
end

%% Get reliability for global metrics
var_global = global_reliability(std_global,nanstd(true_global,0,1));

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
[p,tbl,stats1] = friedman(var_global_80(~isnan(var_global_80(:,1)),:),1,'off');
fprintf('Friedman test for global metrics: p = %1.1e, chi-squared = %1.1f\n',...
    p, tbl{2,5});

% perform a post-hoc Dunn's test
%c = multcompare(stats1,'CType','dunn-sidak');

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

% Compare variability when 80% retained for nodal metrics
[p,tbl,stats1] = friedman(var_nodal_80(~isnan(var_nodal_80(:,1)),:),1,'off');
fprintf('Friedman test for nodal metrics: p = %1.1e, chi-squared = %1.1f\n',...
    p, tbl{2,5});

%c = multcompare(stats1,'CType','dunn-sidak');



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
    %{
    legend('Control centrality','Regional control centrality',...
        'Node strength','Betweenness centrality',...
        'Eigenvector centrality','Clustering coefficient',...
        'location','southeast');
    %}
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
    figure
    set(gcf,'Position',[174 207 1300 850])
    pos_f = get(gcf,'Position');
    [ha, pos] = tight_subplot(3, 2, [0.15 0.04], [-0.05 0.04],[0.05 0.02]);
    set(ha(3),'Position',[pos{3}(1) pos{3}(2)+0.02 ...
        pos{4}(1) + pos{4}(3) - pos{3}(1) pos{3}(4)]);
    set(ha(5),'Position',[pos{5}(1) pos{5}(2)+0.12 ...
        pos{6}(1) + pos{6}(3) - pos{5}(1) pos{5}(4)]);
    delete((ha(4)))
    delete((ha(6)))
    
    % Variability
    axes(ha(1))
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
    xlabel('Percent nodes retained');
    ylabel({'Reliability'})
   % title('Reliability by subsample size','Position',[0.1 0.1 0.1 0.1]);
    set(gca,'Fontsize',20);
    
   % axes(ha(3))
   l1 = legend(nd,{'Control centrality',...
        'Node strength','Betweenness centrality',...
        'Eigenvector centrality','Clustering coefficient'},'Location',...
        'southeast');
    legend boxoff
    
    gl = zeros(size(avg_ag_global,1),1);
    axes(ha(2))
    for j = 1:size(avg_ag_global,1)
         gl(j) = scatter(ef,avg_var_global(j,:),200,'filled');
         hold on
    end
    %{
    legend('Synchronizability','Global efficiency','Transitivity',...
        'location','northeastoutside');
    %}
    xlabel('Percent nodes retained');
   % ylabel({'Reliability'})
    %title('Variability by subsample size');
    set(gca,'Fontsize',20);
    
   % axes(ha(6))
   l2 = legend(gl,{'Synchronizability','Global efficiency','Transitivity'},'Location',...
        'southeast');
    legend boxoff
    annotation('textbox',[0.40 0.91 0.1 0.1],'String',...
        'Reliability by subsample size','FontSize',25,'linestyle','none',...
        'fontweight','bold');
    
    for i = 1:length(l1.EntryContainer.NodeChildren)
        l1.EntryContainer.NodeChildren(i).Icon.Transform.Children.Children.Size = 14;
    end
    
    % Necessary voodoo
    pause(1)
    
    for i = 1:length(l2.EntryContainer.NodeChildren)
        l2.EntryContainer.NodeChildren(i).Icon.Transform.Children.Children.Size = 14;
    end
    

    %% Plot 80% for all patients

    
    cols = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
        0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330];

    % Nodal metrics
    axes(ha(3))
    count = 0;
    for i = 1:size(var_nodal_80,1)
        if isempty(stats(i).name) == 1, continue; end
        count = count + 1;
        for j = 1:size(var_nodal_80,2)
            scatter(count,var_nodal_80(i,j),200,cols(j,:),'filled');
            hold on
        end
    end
   l1 =  legend('Control centrality',...
        'Node strength','Betweenness centrality',...
        'Eigenvector centrality','Clustering coefficient','location',...
        'southeast');
    legend boxoff
    xticks(1:length(names))
    %xlabel('Which patient')
    ylabel('Reliability')
    set(gca,'ylim',[0 1])
    set(gca,'fontsize',20)
    title('Reliability of metric when 20% of network removed','fontsize',25);
    %xticklabels(name_nums)


    % Global metrics
    axes(ha(5))
    count = 0;
    for i = 1:size(var_global_80,1)
        if isempty(stats(i).name) == 1, continue; end
        count = count + 1;
        for j = 1:size(var_global_80,2)
            scatter(count,var_global_80(i,j),200,cols(j,:),'filled');
            hold on
        end
    end
   l2 =  legend('Synchronizability','Global efficiency','Transitivity','location',...
        'southeast');
    legend boxoff
    xticks(1:length(names))
    xlabel('Which patient')
    ylabel('Reliability')
    set(gca,'ylim',[0.7 1])
    set(gca,'fontsize',20)
    
    for i = 1:length(l1.EntryContainer.NodeChildren)
        l1.EntryContainer.NodeChildren(i).Icon.Transform.Children.Children.Size = 14;
    end
    
    % Necessary voodoo
    pause(1)
    
    for i = 1:length(l2.EntryContainer.NodeChildren)
        l2.EntryContainer.NodeChildren(i).Icon.Transform.Children.Children.Size = 14;
    end
    
    annotation('textbox',[0.02 0.86 0.1 0.1],'String',...
        'A','FontSize',35,'linestyle','none');
    annotation('textbox',[0.5 0.86 0.1 0.1],'String',...
        'B','FontSize',35,'linestyle','none');
    annotation('textbox',[0.02 0.49 0.1 0.1],'String',...
        'C','FontSize',35,'linestyle','none');
    annotation('textbox',[0.02 0.2 0.1 0.1],'String',...
        'D','FontSize',35,'linestyle','none');
    
    pause
    print(gcf,[outFolder,'all_fig2_',freq,contig_text,sec_text],'-depsc');
    close(gcf)

end

%xticklabels(name_nums)


end