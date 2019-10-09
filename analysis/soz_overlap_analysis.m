function soz_overlap_analysis(soz_overlap)

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
for freq_idx = 1%1:length(all_freq)
freq = all_freq{freq_idx};

for sec_idx = 1%1:length(all_sec)

    
sec_text = all_sec{sec_idx};

% loop over metrics
all_p = [];
all_z_pts = {};
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
            if isfield(soz_overlap(i).(freq),contig_text) == 0, continue; end

            if isfield(soz_overlap(i).(freq).(contig_text),sec_text) == 0, continue; end

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
        
        
        if 0
        %% Rank sum
        % Comparing single rho from when we only remove soz to all rhos
        % when we randomly remove something that is not the soz
        
        % The reason this is bad is because the group of
        % not_soz_rm_measures are not independent. For instance, there is
        % one patient where the soz electrodes comprise more than half of
        % the electrodes. And so I just take the remainder of the
        % electrodes to remove for not_soz_rm_measure. And I take that same
        % set 1,000 times. And so even a teeny tiny difference will become
        % significant because rank sum assumes that each of the 1,000
        % measures is significant
        
        if exist('soz_rm_measure','var') == 0, continue; end
        [p,h,stats] = ranksum(soz_rm_measure,not_soz_rm_measure);
        
        % Get the z-score
        z = stats.zval;
        
        elseif 1
        %% Difference in means
        % I think this is the only way to directly compare
        if exist('soz_rm_measure','var') == 0, continue; end
        if exist('not_soz_rm_measure','var') == 0, continue; end
        z = (mean(not_soz_rm_measure)-soz_rm_measure);
            
        end
        
        soz_test(metric).name = metrics{metric};
        soz_test(metric).z = [soz_test(metric).z;z];
        soz_test(metric).rm_soz = [soz_test(metric).rm_soz;soz_rm_measure];
        soz_test(metric).rm_not_soz{i} = not_soz_rm_measure;

    end
    
    if isempty(soz_test(metric).z) == 1, continue; end
    
    
    %% T test on the z scores
    [~,p,~,stats] = ttest(soz_test(metric).z);
    soz_test(metric).stats.p = p;
    soz_test(metric).stats.t = stats.tstat;
    soz_test(metric).stats.df = stats.df;
    
    
    % Add to array for plotting
    all_p = [all_p;p];
    all_z_pts{metric} = soz_test(metric).z;

end
 
%% Figure
if isempty(all_p) == 1, continue; end
stars = cell(length(metrics),1);
for i = 1:length(stars)
    if all_p(i) < 0.001/length(metrics)
        stars{i} = '***';
    elseif all_p(i) < 0.01/length(metrics)
        stars{i} = '**';
    elseif all_p(i) < 0.05/length(metrics)
        stars{i} = '*';
    else
        stars{i} = '';
    end
end
cols = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
    0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;...
    0.6350 0.0780 0.1840;0.75 0.5 0.5];


figure
set(gcf,'Position',[174 207 1300 350])
for i = 1:length(all_z_pts)
    scatter(i*ones(size(all_z_pts{i},1),1)+0.05*randn(size(all_z_pts{i},1),1)...
        ,all_z_pts{i},...
       100,'MarkerEdgeColor',cols(i,:),'MarkerFaceColor',cols(i,:))
    hold on
    %{
    scatter(i,all_rho(i),300,'filled','d','MarkerEdgeColor',[0 0.4470 0.7410],...
        'MarkerFaceColor',[0 0.4470 0.7410]);
    %}
    plot([i-0.3,i + 0.3],[mean(all_z_pts{i}),mean(all_z_pts{i})],'color',cols(i,:),'linewidth',3);

    xticks(1:length(all_z_pts))
    xticklabels((metrics))
    xlim([0.7 length(all_z_pts) + 0.3])
    title({'Difference in agreement when ignored electrodes','are seizure onset zone vs non-seizure onset zone'})
    ylabel('Difference in agreement');
    set(gca,'fontsize',20)
    fix_xticklabels(gca,0.1,{'FontSize',20});
end
plot(get(gca,'xlim'),[0 0],'k--','linewidth',2);
for i = 1:length(stars)
    text(i + 0.15, 0.82,stars{i},'fontsize',50);
end


end

end
end

