function compare_soz_resec(soz)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'basic_metrics/'];

contig_text = 'contiguous';
sec_text = 'sec_neg5';
freq = 'high_gamma';


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
dist_names = {'Distance from SOZ','Distance from resection zone',...
    'Overlap with SOZ','Overlap with resection zone',...
    'Average participation coefficient of ignored electrodes',...
    'Average betweenness centrality of ignored electrodes'};
global_metric = [0 0 0 0 0 1 1 1 0 0 0];
    
metrics_to_plot = 6:8;
dist_to_plot = 2;

for dist = dist_to_plot
    figure
    set(gcf,'Position',[8 138 1400 600]);
    [ha, pos] = tight_subplot(min(2,length(metrics_to_plot)/3),...
        min(3,length(metrics_to_plot)), [0.11 0.04], [0.1 0.06],[0.08 0.02]);
    %delete(ha(6));
    count = 0;
    
    for metric = metrics_to_plot
        count = count + 1;
        axes(ha(count))
        
        z = [];
        measure_all = [];
        dist_measure_all = [];
        pt_all = [];
        for i = 1:length(soz)
            if isempty(soz(i).(freq).(contig_text)) == 1, continue; end
            base = soz(i).(freq).(contig_text).(sec_text);
            measure = base.(metrics{metric})';
            
            %% If it's a global metric, take absolute value
            if global_metric(metric) == 1
                measure = abs(measure);
            end
            
            dist_measure = base.(dists{dist})';
            rho =  corr(measure,dist_measure);
            if rho == 1, continue; end
                
            measure_all = [measure_all;measure];
            dist_measure_all = [dist_measure_all;dist_measure];
            pt_all = [pt_all;ones(size(measure,1),1)*i];
            
            z = [z;atanh(rho)];
            
        end
        

        %% Fit a linear mixed effects model
        tbl = table(measure_all,dist_measure_all,pt_all,'VariableNames',...
            {'Metric','Distance','Pt'});
        tbl.Pt = categorical(tbl.Pt);
        lme = fitlme(tbl,'Metric~Distance + (1|Pt)');
        dist_coeff =  lme.Coefficients{2,2};
        dist_p = lme.Coefficients{2,6};
        [~,p] = ttest(z);
        
        
        %% New rho
        new_rho = tanh(nanmean(z));
        
        
        fprintf('For %s and %s, %d had positive correlation and %d had negative.\n\n',...
            metrics{metric},dists{dist},sum(z>0),sum(z<0));
        
        scatter(dist_measure_all,measure_all,50,'DisplayName',sprintf('rho = %1.2f, p = %1.2e',...
            new_rho,p));
        if metric <4
            l = legend('location','southeast');
        else
            l = legend('location','northeast');
        end
        set(l,'Interpreter', 'none');
        l.FontSize = 25;
        
        if count ==1 || count == 4
            if global_metric(metric) == 0
                ylabel({'Spearman rank coefficient','with true metric'});
            elseif global_metric(metric) == 1
                ylabel({'Relative difference','from true metric'});
            end
        end
        if count == 2
            xlabel(sprintf('%s',dist_names{dist}));
        end
        title(sprintf('%s',metric_names{metric}),'Interpreter', 'none');
        set(gca,'fontsize',20)
        print(gcf,[outFolder,dists{dist},'_',freq,contig_text,sec_text],'-depsc');
        %}
        
    end
end
    
   


end