function compare_sz(pt,stats,stats2)

%{
This function looks at the correlation in network metrics from seizure to
seizure and compares this correlation to the correlation from permutation
to permutation
%}

%% Parameters
all_contig = {'random','contiguous'}; % look at random or contiguous removal
all_freq = {'high_gamma','beta'}; % which frequency coherence
all_sec = {'sec_neg10','sec_neg5','sec_0','sec_5','sec_10'}; % which times relative to EEC

%% Locations

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'compare_metrics/'];


%% Initialize arrays
nodal_metrics = {'cc','ns','bc','ec','clust'};
global_metrics = {'sync','eff','trans'};
ef = [20 40 60 80 100];

np = length(stats);

for contig_idx = 1%1:length(all_contig)
for freq_idx = 1%1:length(all_freq)
for sec_idx = 3%1:length(all_sec)
    
    % Get appropriate contig vs random, frequency, time
    contig_text = all_contig{contig_idx};
    sec_text = all_sec{sec_idx};
    freq = all_freq{freq_idx};
    
    for i = 1:length(stats)
    
        if isempty(stats2(i).name) == 1
            if doPlots == 0
                names = [names;nan]; 
            end
            continue; 
        end
    
    
    end
    
    
end
end
end


end