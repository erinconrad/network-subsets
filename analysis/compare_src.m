function compare_src(stats)

%% Parameters
which_freq = 1;
which_contig = 1;
which_sec = 3;

all_contig = {'random','contiguous'}; % look at random or contiguous removal
all_freq = {'high_gamma','beta'}; % which frequency coherence
all_sec = {'sec_neg10','sec_neg5','sec_0','sec_5','sec_10'}; % which times relative to EEC

contig_text = all_contig{which_contig};
sec_text = all_sec{which_sec};
freq = all_freq{which_freq};

nodal_metrics = {'cc','ns','bc','ec','clust'};
nodal_metrics_text = {'control centrality','node strength','betweenness centrality',...
    'eigenvector centrality','clustering coefficient'};
global_metrics = {'sync','eff','trans'};
ef = [20 40 60 80 100];

%% Locations
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'compare_metrics/'];
addpath(genpath(scriptFolder))


%% Initialize arrays of SRC means and variances
src_means = nan(33,5,5); % 33 patients, 5 nodal metrics, 5 removal fractions
src_vars = nan(33,5,5);
names = {};


% Loop through patients
for i = 1:length(stats)
    if isempty(stats(i).name) == 1
        continue; 
    end
    
    %% Extract just numbers from name (for plotting)
    names = [names;stats(i).name];
    [num_idx_s] = regexp(stats(i).name,'\d+');
    %name_nums = [name_nums;stats(i).name(num_idx_s:end)];
    
    if isfield(stats(i).(freq).(contig_text),sec_text) == 0, continue; end
    base = stats(i).(freq).(contig_text).(sec_text);
    
    for m = 1:length(nodal_metrics)
        metric = nodal_metrics{m};
        rho_mean = base.(metric).rho_mean;
        rho_var = base.(metric).rho_var;
        
        % Fill arrays with src mean and var
        src_means(i,m,:) = rho_mean;
        src_vars(i,m,:) = rho_var;
        
    end
    
end

%% Descriptive stats
% I again need to average rhos, but I won't do the Fisher transformation
% here because it is a mess.
src_mean_avg = squeeze(nanmean(src_means,1));
src_mean_std = squeeze(nanstd(src_means,0,1));

%% Statistics to compare nodal metrics
% Do 20% removal
src_means_80 = src_means(:,:,4);
[p,tbl,stats1] = friedman(src_means_80(~isnan(src_means_80(:,1)),:),1,'off');
fprintf('Friedman test for nodal metrics: p = %1.1e, chi-squared = %1.1f, dof = %d\n',...
    p, tbl{2,5},tbl{3,3});

% perform a post-hoc Dunn's test
[c,~,~,~,t] = multcompare_erin(stats1,'CType','dunn-sidak','Display','off');
for i = 1:size(c,1)
    fprintf('Dunn''s test comparing %s and %s: t = %1.2f, p = %1.3f\n\n',...
        nodal_metrics_text{c(i,1)},nodal_metrics_text{c(i,2)},t(i),c(i,6));   
end

%% Prep text
text = cell(5,5);
for f = 1:5
    for m = 1:5
        text{m,f} = sprintf('%1.2f +/- %1.2f',src_mean_avg(m,f),src_mean_std(m,f));
    end
end

table((text(:,1)),(text(:,2)),(text(:,3)),(text(:,4)),'rownames',nodal_metrics)

table((text(:,1)),(text(:,2)),(text(:,3)),(text(:,4)))

end