function compare_soz_resec(soz)

%% Parameters
sec_text = 'sec_neg5';
freq = 'high_gamma';

ex_pts = [1,7]; % example patients to plot correlations for
metrics_to_plot = 1:8;
dist_to_plot = 2;

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
outFolder = [resultsFolder,'basic_metrics/soz/'];
mkdir(outFolder)

measure_1 = [];
measure_2 = [];
dist_1 = [];
dist_2 = [];


for dist = dist_to_plot
    for metric = metrics_to_plot

        z = [];

        for i = 1:length(soz)
            if isempty(soz(i).(freq)) == 1, continue; end
            base = soz(i).(freq).(contig_text).(sec_text);
            
            % Get the agreement metric. For global measures this is the
            % relative difference (which could be positive or negative) and
            % for nodal measures this is the SRC.
            measure = base.(metrics{metric})';
            
            %% If it's a global metric, take absolute value
            % This is because I am just interested in the absolute
            % difference from the true value
            if global_metric(metric) == 1
                measure = abs(measure);
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
        [~,p,stats_t] = ttest(z);
        
        %% Plot figure for this
        f = figure;
        set(f,'Position',[100 100 1400 400]);
        [ha,~] = tight_subplot(1,3,[0.03 0.07],[0.23 0.1],[0.07 0.02]);
        
        % Example patient 1
        axes(ha(1))
        scatter(dist_1,measure_1,100);
        hold on
        [b,~,~,~,stats] = regress(measure_1,[ones(length(dist_1),1),dist_1]);
        r2 = stats(1);
        plot(dist_1,b(1) + b(2)*dist_1)
        
        title(sprintf('Pt 1 - %s',metric_names{metric}))
        xlabel({'Distance of ignored electrodes',dist_names{dist}})
        if global_metric(metric) == 1
            ylabel('Absolute relative difference')
        else
            ylabel('Spearman Rank correlation')
        end
        set(gca,'fontsize',20);
        yl = get(gca,'ylim');
        text(10,yl(1)+0.9*(yl(2)-yl(1)),sprintf('r^2 = %1.2f',r2),'fontsize',20);
        
        % Example patient 2
        axes(ha(2))
        scatter(dist_2,measure_2,100);
        hold on
        [b,~,~,~,stats] = regress(measure_2,[ones(length(dist_2),1),dist_2]);
        r2=stats(1);
        plot(dist_2,b(1) + b(2)*dist_2)
        
        title(sprintf('Pt 2 - %s',metric_names{metric}))
        xlabel({'Distance of ignored electrodes',dist_names{dist}})
        if global_metric(metric) == 1
            ylabel('Absolute relative difference')
        else
            ylabel('Spearman Rank correlation')
        end
        set(gca,'fontsize',20);
        yl = get(gca,'ylim');
        text(10,yl(1)+0.9*(yl(2)-yl(1)),sprintf('r^2 = %1.2f',r2),'fontsize',20);
        
        % All z's
        axes(ha(3))
        scatter(1:length(z),z,200,'filled')
        hold on
        
        xlabel('Which patient')
        ylabel('Fisher''s transformed rho');
        title('Distance leverage association');
        set(gca,'fontsize',20);
        xl = get(gca,'xlim');
        yl = get(gca,'ylim');
        plot(xl,[nanmean(z), nanmean(z)],'linewidth',2);
        text(2,yl(1)+0.9*(yl(2)-yl(1)),sprintf('p = %1.3f',p),'fontsize',20);
        
        pause
        print(gcf,[outFolder,metrics{metric},'_',freq,'_',sec_text,'_',dists{dist}],'-depsc');
        close(gcf)
    end
end
    
   


end