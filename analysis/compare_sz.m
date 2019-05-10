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
outFolder = [resultsFolder,'compare_metrics/'];


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
sz_ag_global = nan(npatients,length(nodal_metrics));

perm_ag_nodal = nan(npatients,length(global_metrics),nf);
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


end