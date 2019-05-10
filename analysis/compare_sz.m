function compare_sz(pt,stats)

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

sz_ag_nodal = nan(np,nf,length(nodal_metrics));
sz_ag_global = nan(np,nf,length(nodal_metrics));

perm_ag_nodal = nan(np,nf,length(nodal_metrics));
perm_ag_global = nan(np,nf,length(nodal_metrics));

for whichPt = npatients

    if isempty(stats(whichPt).name) == 1
        if doPlots == 0
            names = [names;nan]; 
        end
        continue; 
    end
    
    for m = 1:length(all_metrics)
        curr_metric = all_metrics(all_metrics{m});
        base = stats(whichPt).(curr_metric);
        
        sz1 = base.sz(1);
        sz2 = base.sz(2);
        
        if ismember(curr_metric,nodal_metrics) == 1
            
            % rho is agreement measure
            all_rhos_perm = nan(2,nf,np);
            for f = 1:nf
                for p = 1:np
                    all_rhos_perm(1,f,p) = corr(squeeze(sz1.perm(f,p,:)),...
                        sz1.true,'Type','Spearman');
                end
            end
            
        elseif ismember(curr_metric,global_metrics) == 1
            
        end
    end


end
    



end