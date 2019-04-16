function compare_soz_resec(soz,pt)

%{
This function takes patient-level data on agreement between original metric
and resampled metric as a function of distance of the ignored electrodes
from the resection zone and calculates summary statistics and does plots
%}


%% Parameters
doPlot = 0;
all_freq = {'high_gamma','beta'};
all_sec = {'sec_neg10','sec_neg5','sec_0','sec_5','sec_10'};%fieldnames(soz(1).high_gamma.contiguous);
contig_text = 'contiguous'; % Should only do contiguous for this

ex_pts = [1,8]; % example patients to plot correlations for
metrics_to_plot = 1:8; % probably shouldn't change (all nodal and global metrics)
dist_to_plot = 2; % probably shouldn't change (this is distance from resection zone)

metrics = {'rho_cc','rho_ns','rho_bc','rho_ec','rho_clust',...
    'sync','eff','trans',...
    'rho_cc_resec','rho_bc_resec','rho_ns_resec'};
dists = {'dist_soz','dist_resec','overlap_soz','overlap_resec','par_removed',...
    'bc_removed'};
metric_names = {'Control centrality','Node strength','Betweenness centrality',...
    'Eigenvector centrality','Clustering coefficient'...
    'Synchronizability','Global efficiency','Transitivity',...
    'Resection control centrality','Resection betweenness centrality',...
    'Resection node strength'};
dist_names = {'from seizure onset zone','from resection zone',...
    'Overlap with SOZ','Overlap with resection zone',...
    'Average participation coefficient of ignored electrodes',...
    'Average betweenness centrality of ignored electrodes'};
global_metric = [0 0 0 0 0 1 1 1 0 0 0];

%% Locations
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'agreement_dist/'];

%% Initialize cross time, freq, etc. comparisons
t_all = nan(length(metrics_to_plot),length(all_sec),...
    length(all_freq));
p_all = nan(length(metrics_to_plot),length(all_sec),...
    length(all_freq));
t_text = cell(length(metrics_to_plot),length(all_sec),...
    length(all_freq));


%{
% Should only do contiguous for this

sec_text = 'sec_0';
freq = 'high_gamma';
    %}


for freq_idx = 1%1:length(all_freq)
for sec_idx = 3%1:length(all_sec)
    

sec_text = all_sec{sec_idx};
freq = all_freq{freq_idx};
measure_1 = [];
measure_2 = [];
dist_1 = [];
dist_2 = [];
name_1 = pt(ex_pts(1)).name;
name_2 = pt(ex_pts(2)).name;
all_p = [];
all_t = [];
all_rho = [];
all_z = [];
all_rho_pts = {};

for dist = dist_to_plot
    for metric = metrics_to_plot

        z = [];
        rho_pts = [];

        for i = 1:length(soz)
            if isempty(soz(i).(freq)) == 1, continue; end
            if isfield(soz(i).(freq).(contig_text),sec_text) == 0, continue; end
            base = soz(i).(freq).(contig_text).(sec_text);
            
            % Get the agreement metric. For global measures this is the
            % relative difference (which could be positive or negative) and
            % for nodal measures this is the SRC.
            measure = base.(metrics{metric})';
            
            %% If it's a global metric, take absolute value and make negative
            % This is because I am just interested in the absolute
            % difference from the true value; negative to make it go in the
            % same direction as nodal measure (higher = more agreement)
            if global_metric(metric) == 1
                measure = -abs(measure);
            end
            
            % Get the measure of distance (usually doing distance from
            % resection zone)
            dist_measure = base.(dists{dist})';
            
            % Correlate the agreement metric with the distance metric
            rho =  corr(measure,dist_measure,'Type','Spearman');
            
            % I am not sure why this would be the case
            if rho == 1, error('what\n'); end
 
            % Aggregate the transformed rho's for each patient
            z = [z;atanh(rho)];
            rho_pts = [rho_pts;rho];
            
            
            if i == ex_pts(1)
                measure_1 = measure;
                dist_1 = dist_measure;
                
            elseif i == ex_pts(2)
                measure_2 = measure;
                dist_2 = dist_measure;
            end
            if 1 == 0
            figure
            scatter(dist_measure,measure)
            title(metrics{metric});
            pause
            close(gcf)
            end
        end
        
        
        %% T test to see if the aggregated transformed rho's diff from zero
        
        
        all_rho = [all_rho;average_rho(rho_pts,1)];
        
        % two sided one sample t test to see if the transformed rho's,
        % aggregated across patients, are different from zero
        [~,p,~,stats] = ttest(z); 
        all_z = [all_z;nanmean(z)];
        all_p = [all_p;p];
        all_t = [all_t;stats.tstat];
        all_rho_pts = [all_rho_pts;rho_pts];
        
        %% Plot figure for this
        %{
        f = figure;
        set(f,'Position',[100 100 1450 400]);
        [ha,~] = tight_subplot(1,3,[0.03 0.07],[0.23 0.1],[0.05 0.03]);
        
        % Example patient 1
        axes(ha(1))
        scatter(dist_1,measure_1,100,'linewidth',2);
        hold on
        [b,~,~,~,stats] = regress(measure_1,[ones(length(dist_1),1),dist_1]);
        r2 = stats(1);
        plot(dist_1,b(1) + b(2)*dist_1)
        
        title(sprintf('%s - %s',metric_names{metric},name_1))
        if dist == 2
            xlabel({'Distance of ignored electrodes',dist_names{dist}})
        elseif dist == 6
            xlabel({'Average betweenness centrality', 'of ignored electrodes'})
        end
        if global_metric(metric) == 1
            ylabel('Absolute relative difference')
        else
            ylabel('Spearman''s rank correlation')
        end
        set(gca,'fontsize',20);
        yl = get(gca,'ylim');
        text(15,yl(1)+0.1*(yl(2)-yl(1)),sprintf('r^2 = %1.2f',r2),'fontsize',20,...
            'fontweight','bold');
        
        % Example patient 2
        axes(ha(2))
        scatter(dist_2,measure_2,100,'linewidth',2);
        hold on
        [b,~,~,~,stats] = regress(measure_2,[ones(length(dist_2),1),dist_2]);
        r2=stats(1);
        plot(dist_2,b(1) + b(2)*dist_2)
        
        title(sprintf('%s - %s',metric_names{metric},name_2))
        if dist == 2
            xlabel({'Distance of ignored electrodes',dist_names{dist}})
        elseif dist == 6
            xlabel({'Average betweenness centrality', 'of ignored electrodes'})
        end
        if global_metric(metric) == 1
            ylabel('Absolute relative difference')
        else
            ylabel('Spearman''s rank correlation')
        end
        set(gca,'fontsize',20);
        yl = get(gca,'ylim');
        text(15,yl(1)+0.1*(yl(2)-yl(1)),sprintf('r^2 = %1.2f',r2),'fontsize',20,...
            'fontweight','bold');
        
        % All z's
        axes(ha(3))
        scatter(1:length(z),z,200,'linewidth',2)
        hold on
        
        xlabel('Which patient')
        ylabel('Fisher''s transformed rho');
        if dist == 2
            title('Distance association - all patients');
        elseif dist == 6
            title('Betweenness centrality association - all');
        end
        set(gca,'fontsize',20);
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        plot(xl,[nanmean(z), nanmean(z)],'linewidth',2);
        if p < 0.001
            text(2,yl(1)+0.9*(yl(2)-yl(1)),'p < 0.001','fontsize',20,...
            'fontweight','bold');
        else
            text(2,yl(1)+0.9*(yl(2)-yl(1)),sprintf('p = %1.3f',p),'fontsize',20,...
            'fontweight','bold');
        end
        
        pause
        print(gcf,[outFolder,metrics{metric},'_',freq,'_',sec_text,'_',dists{dist}],'-depsc');
        close(gcf)
        %}
    end
end

t_all(:,sec_idx,freq_idx) = all_t;
p_all(:,sec_idx,freq_idx) = all_p;
for i = 1:size(all_t)
    if all_p(i) < 0.001/length(metrics_to_plot)
        extra = '***';
    elseif all_p(i) < 0.01/length(metrics_to_plot)
        extra = '**';
    elseif all_p(i) < 0.05/length(metrics_to_plot)
        extra = '*';
    else
        extra = '';
    end
    t_text{i,sec_idx,freq_idx} = sprintf('%1.2f%s',all_t(i),extra);
end

end
end
%{
%% Table with different times, high gamma
%
table(char(t_text(:,1,1)),char(t_text(:,2,1)),char(t_text(:,3,1)),...
    char(t_text(:,4,1)),char(t_text(:,5,1)),'VariableNames',all_sec,...
    'RowNames',metrics(1:8))
%}

% For sec_0, high_gamma, plot t scores and p values
%{
table(t_text(:,3,1),p_all(:,3,1),'RowNames',metrics(1:8))
%}

%{
%% EEC, beta
squeeze(t_all(:,3,2))
squeeze(p_all(:,3,2))

%% EEC, high gamma (for sz 2)
squeeze(t_all(:,3,1))
squeeze(p_all(:,3,1))
%}


%% Set high_gamma and compare times

if doPlot == 1

    
    sig_locs_x = find(all_p < 0.05);
    sig_locs_y = 0.9*ones(length(sig_locs_x),1);
    stars = cell(length(metrics_to_plot),1);
    for i = 1:length(stars)
        if all_p(i) < 0.001/length(metrics_to_plot)
            stars{i} = '***';
        elseif all_p(i) < 0.01/length(metrics_to_plot)
            stars{i} = '**';
        elseif all_p(i) < 0.05/length(metrics_to_plot)
            stars{i} = '*';
        else
            stars{i} = '';
        end
    end
    cols = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
        0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;...
        0.6350 0.0780 0.1840;0.75 0.5 0.5];

    
    figure
    set(gcf,'Position',[174 207 1300 350])
    for i = 1:length(all_rho_pts)
        scatter(i*ones(size(all_rho_pts{i},1),1)+0.05*randn(size(all_rho_pts{i},1),1)...
            ,all_rho_pts{i},...
           100,'MarkerEdgeColor',cols(i,:),'MarkerFaceColor',cols(i,:))
        hold on
        %{
        scatter(i,all_rho(i),300,'filled','d','MarkerEdgeColor',[0 0.4470 0.7410],...
            'MarkerFaceColor',[0 0.4470 0.7410]);
        %}
        plot([i-0.3,i + 0.3],[all_rho(i),all_rho(i)],'color',cols(i,:),'linewidth',3);

        xticks(1:length(all_rho_pts))
        xticklabels(metric_names(metrics_to_plot))
        xlim([0.7 length(all_rho) + 0.3])
        title({'Association between metric accuracy and','distance of ignored electrodes from resection zone'})
        ylabel('Distance-agreement correlation');
        set(gca,'fontsize',20)
        fix_xticklabels(gca,0.1,{'FontSize',20});
    end
    plot(get(gca,'xlim'),[0 0],'k--','linewidth',2);
    for i = 1:length(stars)
        text(i + 0.15, 0.82,stars{i},'fontsize',50);
    end
    pause
    print(gcf,[outFolder,'all_',freq,'_',sec_text,'_',dists{dist}],'-depsc');
    close(gcf)

    % convert all_rho_pts to matrix
    M = cell2mat(all_rho_pts');

    %% Box plot
    figure
    set(gcf,'Position',[174 207 1300 350])
    boxplot(M)
    hold on
    xticks(1:length(all_rho_pts))
    xticklabels(metric_names(metrics_to_plot))
    xlim([0.7 length(all_rho) + 0.3])
    title({'Association between metric accuracy and','distance of ignored electrodes from resection zone'})
    ylabel('Distance-agreement correlation');
    set(gca,'fontsize',20)
    fix_xticklabels(gca,0.1,{'FontSize',20});
    plot(get(gca,'xlim'),[0 0],'color',[0.8500 0.3250 0.0980],'linewidth',2);
    scatter(sig_locs_x + 0.15, sig_locs_y + 0.01,300,'*','k');
    pause
    print(gcf,[outFolder,'box_plot_',freq,'_',sec_text,'_',dists{dist}],'-depsc');
    close(gcf)
%}

     % convert all_rho_pts to matrix
    M = cell2mat(all_rho_pts');
    %% Violin plot
    figure
    set(gcf,'Position',[174 207 1300 350])
    violin(M,1)
    hold on
    xticks(1:length(all_rho_pts))
    xticklabels(metric_names(metrics_to_plot))
    xlim([0.7 length(all_rho) + 0.3])
    title({'Association between metric accuracy and','distance of ignored electrodes from resection zone'})
    ylabel('Distance-agreement correlation');
    set(gca,'fontsize',20)
    fix_xticklabels(gca,0.1,{'FontSize',20});
    plot(get(gca,'xlim'),[0 0],'k--','linewidth',2);
    for i = 1:length(stars)
        text(i + 0.15, 0.9 + 0.08,stars{i},'fontsize',50);
    end
    ax = gca;
    outerpos = ax.OuterPosition;
    ti = ax.TightInset; 
    left = outerpos(1) + ti(1);
    bottom = outerpos(2) + ti(2);
    ax_width = outerpos(3) - ti(1) - ti(3);
    ax_height = outerpos(4) - ti(2) - ti(4);
    ax.Position = [left+0.01 bottom + 0.1 ax_width-0.02 ax_height-0.1];
    pause
    print(gcf,[outFolder,'violin_plot_',freq,'_',sec_text,'_',dists{dist}],'-depsc');
    print(gcf,[outFolder,'violin_plot_',freq,'_',sec_text,'_',dists{dist}],'-dpng');
    close(gcf)

end

table(metrics(metrics_to_plot)',all_rho,all_t,all_p)



%{
figure
set(gcf,'Position',[174 207 1300 350])
stem(1:length(all_rho),all_rho,'filled','markersize',15)
hold on
scatter(sig_locs_x + 0.15, sig_locs_y + 0.01,150,'*','k');
xticklabels(metric_names(metrics_to_plot))
xlim([0.7 length(all_rho) + 0.3])
title({'Association between metric accuracy and','distance of ignored electrodes from resection zone'})
ylabel('Distance-agreement correlation');
set(gca,'fontsize',20)
fix_xticklabels(gca,0.1,{'FontSize',20});
print(gcf,[outFolder,'all_',freq,'_',sec_text,'_',dists{dist}],'-depsc');
%}

end