function stats_CI(stats)

%{
This function takes info about the number and identities of electrodes
forming the nth% CI for highest metric value and calculates summary stats about
them
%}

%% Parameters
all_contig = {'random','contiguous'};
all_freq = {'high_gamma','beta'};
all_sec = {'sec_neg10','sec_neg5','sec_0','sec_5','sec_10'};
nodal_metrics = {'ns','bc','ec','clust','cc','regional_cc'};
global_metrics = {'sync','eff','trans'};

all_global = nan(length(global_metrics),length(all_sec),...
    length(all_freq),length(all_contig));
all_nodal = nan(length(nodal_metrics),length(all_sec),...
    length(all_freq),length(all_contig));

%% Initialize parameters
np = length(stats);

for contig_idx = 1:length(all_contig)
for freq_idx = 1:length(all_freq)
for sec_idx = 1:length(all_sec)
    
contig_text = all_contig{contig_idx};
sec_text = all_sec{sec_idx};
freq = all_freq{freq_idx};

true_nums_nodal = nan(np,length(nodal_metrics));
nums_nodal_95 = nan(np,length(nodal_metrics));
global_width_95 = nan(np,length(global_metrics));
names = cell(np,1);

%% Get true numbers and 95% CI info 
for i = 1:length(stats)
    
    if isempty(stats(i).name) == 1, continue; end
    if isfield(stats(i).(freq).(contig_text),sec_text) == 0, continue; end
    base = stats(i).(freq).(contig_text).(sec_text);
    
    names{i} = stats(i).name;
    
    % Get all nodal metrics other than regional cc
    for j = 1:length(nodal_metrics)-1
        true_nums_nodal(i,j) = 1;
        nums_nodal_95(i,j) = length(base.(nodal_metrics{j}).single_95);
    end
    
    % Get regional cc
    if isnan(base.min_cc_elecs.regional_95) == 0
        true_nums_nodal(i,end) = length(base.regional_cc.true);
        nums_nodal_95(i,end) = length(base.min_cc_elecs.regional_95);
    end
    
    % Get all global metrics
    for j = 1:length(global_metrics)
        t = base.(global_metrics{j}).all;
        w = prctile(t,[2 97.5]);
        global_width_95(i,j) = w(2)-w(1);
    end
    
end

%% Get stats for number and ratio of channels in 95% CI for each nodal metric
% This can be presented as a mean (min-max) for each metric. For regional
% cc, can present both number and ratio.
for j = 1:length(nodal_metrics)
    fprintf(['The median number of electrodes in the 95%% CI for\n%s is '...
        '%1.1f (range %d-%d).\n\n'],nodal_metrics{j},nanmedian(nums_nodal_95(:,j)),...
        min(nums_nodal_95(:,j)),max(nums_nodal_95(:,j)));
end

fprintf(['The median ratio of electrodes in the 95%% CI for\n%s is '...
        '%1.2f (range %1.2f-%1.2f).\n\n'],nodal_metrics{end},...
        nanmedian(nums_nodal_95(:,end)./true_nums_nodal(:,end)),...
        min(nums_nodal_95(:,end)./true_nums_nodal(:,end)),...
        max(nums_nodal_95(:,end)./true_nums_nodal(:,end)));


%% Get stats for 95% CI for each global metric
% This can be presented as mean (min-max) of the width of the 95% CI for
% each metric.

for j = 1:length(global_metrics)
    fprintf(['The median width of the 95%% CI for %s is\n%1.3f '...
        '(range %1.3f-%1.3f).\n\n'],global_metrics{j},nanmedian(global_width_95(:,j)),...
        min(global_width_95(:,j)),max(global_width_95(:,j)));
    
end

all_global(:,sec_idx,freq_idx,contig_idx) = nanmedian(global_width_95,1);
all_nodal(:,sec_idx,freq_idx,contig_idx) = nanmedian(nums_nodal_95,1);

%{
%% Table to probe nodal measures
table(names,nums_nodal_95(:,1),nums_nodal_95(:,2),nums_nodal_95(:,3),...
    nums_nodal_95(:,4),nums_nodal_95(:,5),nums_nodal_95(:,6),'VariableNames',...
    {'names',nodal_metrics{1},nodal_metrics{2},nodal_metrics{3},nodal_metrics{4},...
    nodal_metrics{5},nodal_metrics{6}})

%% Table to probe global measures
table(names,global_width_95(:,1),global_width_95(:,2),global_width_95(:,3))
%}
end
end
end

%% Vary times
all_nodal(:,:,1,1)
all_global(:,:,1,1)

%% Beta
all_nodal(:,3,2,1)
all_global(:,3,2,1)

%% Contig
all_nodal(:,3,1,2)
all_global(:,3,1,2)

%% EEC, high gamma, random (sz 2)
all_nodal(:,3,1,1)
all_global(:,3,1,1)


end