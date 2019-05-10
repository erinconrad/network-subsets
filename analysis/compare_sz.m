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
ef = [20 40 60 80 100];
nf = length(ef);
np = length(stats);
names = [];

sz_ag_nodal = nan(np,nf,length(nodal_metrics));

for i = 1:length(stats)

    if isempty(stats(i).name) == 1
        if doPlots == 0
            names = [names;nan]; 
        end
        continue; 
    end
    
    % Get \


end
    



end