function soz_overlap_analysis(soz_overlap,pt)

%% Parameters
metrics = {'rho_cc','rho_ns','rho_bc','rho_ec','rho_clust',...
    'sync','eff','trans'};
n_metrics = length(metrics);
global_metric = [0 0 0 0 0 1 1 1];
all_freq = {'high_gamma','beta'};
all_sec = {'sec_neg10','sec_neg5','sec_0','sec_5','sec_10'};
all_contig = {'not_soz','soz'};

%% Locations
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'soz_overlap/'];

% Loop over frequencies and times
for freq_idx = 1:length(all_freq)
freq = all_freq{freq_idx};

for sec_idx = 1:length(all_sec)
skip_flag = 0;
    
sec_text = all_sec{sec_idx};

% loop over metrics
for metric = 1:n_metrics
    
    % Initialize some variables
    soz_test(metric).z = [];
    soz_test(metric).rm_soz = [];
    soz_test(metric).rm_not_soz = {};

    % loop over patients
    for i = 1:length(soz_overlap)
        if isfield(soz_overlap(i),freq) == 0, continue; end
        if isempty(soz_overlap(i).(freq)) == 1, continue; end
        
        
        
        for contig_idx = [1,2]
    
            contig_text = all_contig{contig_idx};
        
            if isfield(soz_overlap(i).(freq).(contig_text),sec_text) == 0
                skip_flag = 1;
                continue; 
            end

            % Get base
            base = soz_overlap(i).(freq).(contig_text).(sec_text);

            % Get agreement metric
            measure = base.(metrics{metric})';

            % If it's a global metric, take absolute value and make negative
            % This is because I am just interested in the absolute
            % difference from the true value; negative to make it go in the
            % same direction as nodal measure (higher = more agreement)
            if global_metric(metric) == 1
                measure = -abs(measure);
            end
            
            if contig_idx == 2
                soz_rm_measure = measure;
            elseif contig_idx == 1
                not_soz_rm_measure = measure;
            end
        
        end
        
        if skip_flag == 1, continue; end
        
        %% Rank sum
        % Comparing single rho from when we only remove soz to all rhos
        % when we randomly remove something that is not the soz
        [p,h,stats] = ranksum(soz_rm_measure,not_soz_rm_measure);
        
        % Get the z-score
        z = stats.zval;
        
        soz_test(metric).name = metrics{metric};
        soz_test(metric).z = [soz_test(metric).z;z];
        soz_test(metric).rm_soz = [soz_test(metric).rm_soz;soz_rm_measure];
        soz_test(metric).rm_not_soz{i} = not_soz_rm_measure;

    end

end
 
if isempty(soz_test(1).z) == 0
    error('look\n');
end

end

end
end

