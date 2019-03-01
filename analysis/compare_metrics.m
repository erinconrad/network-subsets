function compare_metrics(stats)


[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'basic_metrics/'];

doPlots = 1;

%% Initialize arrays
nodal_metrics = {'cc','ns','bc'};
global_metrics = {'sync','eff'};
ef = [20 40 60 80 100];

np = length(stats);

names ={};
name_nums = {};
ag_nodal = zeros(np,length(nodal_metrics),length(ef));
var_nodal = zeros(np,length(nodal_metrics),length(ef));
ag_global = zeros(np,length(global_metrics),length(ef));
std_global = zeros(np,length(global_metrics),length(ef));
true_global = zeros(np,length(global_metrics));

for i = 1:length(stats)
    
    %% Extract just numbers from name (for plotting)
    names = [names;stats(i).name];
    [num_idx_s] = regexp(stats(i).name,'\d+');
    name_nums = [name_nums;stats(i).name(num_idx_s:end)];
    
    %% Look at contiguous and -5 seconds
    contig_text = 'random';
    sec_text = 'sec_neg5';
    
    base = stats(i).(contig_text).(sec_text);
    
    %% Get agreement and variability for metrics by resection size
    for j = 1:length(nodal_metrics)
        var_nodal(i,j,:) = base.(nodal_metrics{j}).rel_std;
        ag_nodal(i,j,:) = base.(nodal_metrics{j}).rho_mean';
    end
    
    %% Get agreement and variability for global metrics by resection size
    for j = 1:length(global_metrics)
        std_global(i,j,:) = base.(global_metrics{j}).std';
        ag_global(i,j,:) = mean(base.(global_metrics{j}).rel_diff,2);
        true_global(i,j,:) = base.(global_metrics{j}).true;
    end
   
end

%% Get var_global (relative std by dividing by std across patients)
var_global = std_global./std(true_global,0,1);

%% Average over patients
avg_ag_nodal = squeeze(average_rho(ag_nodal,1)); % Fisher transform for rho
avg_var_nodal = squeeze(mean(var_nodal,1));
avg_ag_global = squeeze(mean(ag_global,1));
avg_var_global = squeeze(mean(var_global,1));

%% Individual patient, 80% retained
var_nodal_80 = var_nodal(:,:,4);
var_global_80 = var_global(:,:,4);

%% Calculate statistics for variability

% Compare variability when 80% retained for global metrics
[p,~,stats] = signrank(var_global_80(:,1),var_global_80(:,2));
fprintf('Sign rank test for global metrics: p = %1.2e, sign-rank = %d\n',...
    p, stats.signedrank);

% Compare variability when 80% retained for nodal metrics
[p,tbl] = friedman(var_nodal_80,1,'off');
fprintf('Friedman test for global metrics: p = %1.1e, chi-squared = %1.1f\n',...
    p, tbl{2,5});

[p,~,stats] = signrank(var_nodal_80(:,1),var_nodal_80(:,2));
fprintf('Sign rank test for cc vs ns: p = %1.2e, sign-rank = %d\n',...
    p, stats.signedrank);

[p,~,stats] = signrank(var_nodal_80(:,1),var_nodal_80(:,3));
fprintf('Sign rank test for cc vs bc: p = %1.2e, sign-rank = %d\n',...
    p, stats.signedrank);

[p,~,stats] = signrank(var_nodal_80(:,2),var_nodal_80(:,3));
fprintf('Sign rank test for ns vs bc: p = %1.2e, sign-rank = %d\n',...
    p, stats.signedrank);

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
    legend('Control centrality','Node strength','Betweenness centrality',...
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
    legend('Control centrality','Node strength','Betweenness centrality',...
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
    legend('Synchronizability','Global efficiency',...
        'location','southeast');
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
    legend('Synchronizability','Global efficiency',...
        'location','northeast');
    xlabel('Percent nodes retained');
    ylabel({'Relative standard deviation'})
    %title('Variability by subsample size');
    set(gca,'Fontsize',20);
    pause
    print(gcf,[outFolder,'avg_metrics_',contig_text,sec_text],'-depsc');
    close(gcf)

    %% Plot 80% for all patients

    figure
    set(gcf,'Position',[31 173 1400 632]);
    [ha, pos] = tight_subplot(2, 1, [0.1 0.15], [0.1 0.07],[0.06 0.02]);
    cols = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250];

    % Nodal metrics
    axes(ha(1))
    for i = 1:size(var_nodal_80,1)
        for j = 1:size(var_nodal_80,2)
            scatter(i,var_nodal_80(i,j),200,cols(j,:),'filled');
            hold on
        end
    end
    legend('Control centrality','Node strength','Betweenness centrality');
    xticks(1:length(names))
    %xlabel('Which patient')
    ylabel('Relative standard deviation')
    set(gca,'fontsize',20)
    title('Variability of metric when 20% of network removed');
    %xticklabels(name_nums)


    % Global metrics
    axes(ha(2))
    for i = 1:size(var_global_80,1)
        for j = 1:size(var_global_80,2)
            scatter(i,var_global_80(i,j),200,cols(j,:),'filled');
            hold on
        end
    end
    legend('Synchronizability','Global efficiency');
    xticks(1:length(names))
    xlabel('Which patient')
    ylabel('Relative standard deviation')
    set(gca,'fontsize',20)
    pause
    print(gcf,[outFolder,'all_pat_80_',contig_text,sec_text],'-depsc');
    close(gcf)

end

%xticklabels(name_nums)


end