function compare_soz_resec(soz,pt)

%% Parameters
sec_text = 'sec_neg5';
freq = 'high_gamma';

ex_pts = [1,8]; % example patients to plot correlations for
metrics_to_plot = 1:8;%1:8;
dist_to_plot = 2; % probably shouldn't change (this is distance from resection zone)

metrics = {'rho_cc','rho_bc','rho_ns','rho_ec','rho_clust',...
    'eff','sync','trans',...
    'rho_cc_resec','rho_bc_resec','rho_ns_resec'};
dists = {'dist_soz','dist_resec','overlap_soz','overlap_resec','par_removed',...
    'bc_removed'};
metric_names = {'Control centrality','Betweenness centrality','Node strength',...
    'Eigenvector centrality','Clustering coefficient'...
    'Global efficiency','Synchronizability','Transitivity',...
    'Resection control centrality','Resection betweenness centrality',...
    'Resection node strength'};
dist_names = {'from seizure onset zone','from resection zone',...
    'Overlap with SOZ','Overlap with resection zone',...
    'Average participation coefficient of ignored electrodes',...
    'Average betweenness centrality of ignored electrodes'};
global_metric = [0 0 0 0 0 1 1 1 0 0 0];


%% Don't change
% Should only do contiguous for this
contig_text = 'contiguous';

%% Locations

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'agreement_dist/'];

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
            base = soz(i).(freq).(contig_text).(sec_text);
            
            % Get the agreement metric. For global measures this is the
            % relative difference (which could be positive or negative) and
            % for nodal measures this is the SRC.
            measure = base.(metrics{metric})';
            
            %% If it's a global metric, take absolute value
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
            rho =  corr(measure,dist_measure);
            
            % I am not sure why this would be the case
            if rho == 1, continue; end
 
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

sig_locs_x = find(all_p < 0.05);
sig_locs_y = 0.9*ones(length(sig_locs_x),1);

figure
set(gcf,'Position',[174 207 1300 350])
for i = 1:length(all_rho_pts)
    scatter(i*ones(size(all_rho_pts{i},1),1)+0.05*randn(size(all_rho_pts{i},1),1)...
        ,all_rho_pts{i},...
       100,'MarkerEdgeColor',[0 0.4470 0.7410])
    hold on
    scatter(i,all_rho(i),300,'filled','d','MarkerEdgeColor',[0 0.4470 0.7410],...
        'MarkerFaceColor',[0 0.4470 0.7410]);
    
    xticks(1:length(all_rho_pts))
    xticklabels(metric_names(metrics_to_plot))
    xlim([0.7 length(all_rho) + 0.3])
    title({'Association between metric accuracy and','distance of ignored electrodes from resection zone'})
    ylabel('Distance-agreement correlation');
    set(gca,'fontsize',20)
    fix_xticklabels(gca,0.1,{'FontSize',20});
end
plot(get(gca,'xlim'),[0 0],'color',[0.8500 0.3250 0.0980],'linewidth',2);
scatter(sig_locs_x + 0.15, sig_locs_y + 0.01,200,'*','k');
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