function compare_sz(stats)

%{
This function looks at the correlation in network metrics from seizure to
seizure and compares this correlation to the correlation from permutation
to permutation
%}

%% Parameters

%% Locations

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'sz_comp/'];


%% Initialize arrays
nodal_metrics = {'cc','ns','bc','ec','clust'};
global_metrics = {'sync','eff','trans'};
all_metrics = [nodal_metrics,global_metrics];
ef = [20 40 60 80 100];
nf = length(ef);
np = 1000;
npatients = length(stats);
names = [];

sz_ag_nodal = nan(npatients,length(nodal_metrics));
sz_ag_global = nan(npatients,length(global_metrics));

perm_ag_nodal = nan(npatients,length(nodal_metrics),nf);
perm_ag_global = nan(npatients,length(global_metrics),nf);

for whichPt = 1:npatients

    if isempty(stats(whichPt).cc) == 1
        
        continue; 
    end
    
    for m = 1:length(all_metrics)
        curr_metric = all_metrics{m};
        base = stats(whichPt).(curr_metric);
        
        sz1 = base.sz(1);
        sz2 = base.sz(2);
        
        if ismember(curr_metric,nodal_metrics) == 1
            
            % rho is agreement measure - get rho between each perm and true
            all_rhos_perm = nan(2,nf,np);
            for f = 2:nf
                for p = 1:np
                    all_rhos_perm(1,f,p) = corr(squeeze(sz1.perm(f,p,...
                        ~isnan(sz1.perm(f,p,:))&~isnan(sz1.perm(f,1,:)))),...
                        squeeze(sz1.perm(f,1,...
                        ~isnan(sz1.perm(f,p,:))&~isnan(sz1.perm(f,1,:)))),...
                        'Type','Spearman');
                    
                    all_rhos_perm(2,f,p) = corr(squeeze(sz2.perm(f,p,...
                        ~isnan(sz2.perm(f,p,:))&~isnan(sz2.perm(f,1,:)))),...
                        squeeze(sz2.perm(f,1,...
                        ~isnan(sz2.perm(f,p,:))&~isnan(sz2.perm(f,1,:)))),...
                        'Type','Spearman');
                end
                
                
            end
            
            % Get sz rho
            sz_ag_nodal(whichPt,m) = corr(sz1.true,sz2.true,'Type','Spearman');
            
            % Do average of rhos across all perms and szs
            avg_rhos_perm = squeeze(mean(all_rhos_perm,3));
            avg_rhos_perm = squeeze(mean(avg_rhos_perm,1));
            perm_ag_nodal(whichPt,m,:) = avg_rhos_perm;
               
            
        elseif ismember(curr_metric,global_metrics) == 1
            
            % -absolute value of relative diff is agreement measure
            all_aggs_perm = nan(2,nf,np);
            for f = 1:nf
                for p = 1:np
            
                    all_aggs_perm(1,f,p) = -abs((squeeze(sz1.perm(f,p))-...
                        squeeze(sz1.perm(f,1)))/squeeze(sz1.perm(f,1)));
                    
                    all_aggs_perm(2,f,p) = -abs((squeeze(sz2.perm(f,p))-...
                        squeeze(sz2.perm(f,1)))/squeeze(sz2.perm(f,1)));
            
                end
                
               
                
            end
            
             % Get sz agg
            sz_ag_global(whichPt,m-length(nodal_metrics)) = -abs((sz2.true-sz1.true)/sz1.true);
            
            % Do average of ags across perms and szs
            avg_aggs_perm = squeeze(nanmean(nanmean(all_aggs_perm,3),1));
            perm_ag_global(whichPt,m-length(nodal_metrics),:) = avg_aggs_perm;
            
        end
    end


end
    
%% For 20% removal, do a paired ttest for each metric to see which better
cols = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
        0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;...
        0.6350 0.0780 0.1840;0.75 0.5 0.5];
figure
set(gcf,'position',[100 100 800 700])
[ha, pos] = tight_subplot(2, 1, [0.13 0.04], [0.05 0.05],[0.09 0.02]);
axes(ha(1))
xtickpos = [];
xticknames = {};
for j = 1:length(nodal_metrics)
    er = squeeze(nanstd(sz_ag_nodal,0,1));
    errorbar(j,squeeze(nanmean(sz_ag_nodal(:,j),1)),er(j),...
        'o','MarkerSize',15,'MarkerEdgeColor',cols(j,:),...
       'MarkerFaceColor',cols(j,:),'Color',cols(j,:),'linewidth',2);
    hold on
    
    errorbar(j+0.3,squeeze(nanmean(perm_ag_nodal(:,j,4),1)),er(j),...
        'd','MarkerSize',15,'MarkerEdgeColor',cols(j,:),...
       'MarkerFaceColor',cols(j,:),'Color',cols(j,:),'linewidth',2); 
   
   xtickpos = [xtickpos j j+0.3];
   xticknames = [xticknames, sprintf('%s sz',nodal_metrics{j}),...
       sprintf('%s perm',nodal_metrics{j})];
end
xticks(xtickpos);
xticklabels(xticknames);
fix_xticklabels(gca,0.1,{'FontSize',20});
ylabel('Spearman rank correlation')
title('Nodal metric agreement')
set(gca,'fontsize',20)

xtickpos = [];
xticknames = {};
axes(ha(2))
for j = 1:length(global_metrics)
    er = squeeze(nanstd(sz_ag_global,0,1));
    errorbar(j,squeeze(nanmean(sz_ag_global(:,j),1)),er(j),...
        'o','MarkerSize',15,'MarkerEdgeColor',cols(j+5,:),...
       'MarkerFaceColor',cols(j+5,:),'Color',cols(j+5,:),'linewidth',2);
    hold on
    
    errorbar(j+0.3,squeeze(nanmean(perm_ag_global(:,j,4),1)),er(j),...
        'd','MarkerSize',15,'MarkerEdgeColor',cols(j+5,:),...
       'MarkerFaceColor',cols(j+5,:),'Color',cols(j+5,:),'linewidth',2); 
   
   xtickpos = [xtickpos j j+0.3];
   xticknames = [xticknames, sprintf('%s sz',global_metrics{j}),...
       sprintf('%s perm',global_metrics{j})];
end
xticks(xtickpos);
xticklabels(xticknames);
fix_xticklabels(gca,0.1,{'FontSize',20});
ylabel('Negative relative difference')
title('Global metric agreement')
set(gca,'fontsize',20)
print(gcf,[outFolder,'thing'],'-depsc')

end