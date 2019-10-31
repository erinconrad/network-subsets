function soz_overlap_analysis(soz_overlap,pt)


%{
This is the analysis we added for Network Neuroscience. It compares metrics
when we target vs spare the SOZ for removal
%}

%% Parameters
do_plot = 0;
save_plot = 0;
metrics = {'rho_cc','rho_ns','rho_bc','rho_ec','rho_clust',...
    'sync','eff','trans'};
hub_metrics = {'cc','ns','bc','ec','clust'};
metric_names = {'Control centrality','Node strength','Betweenness centrality',...
    'Eigenvector centrality','Clustering coefficient'...
    'Synchronizability','Global efficiency','Transitivity'};

n_metrics = length(metrics);
global_metric = [0 0 0 0 0 1 1 1];
all_freq = {'high_gamma','beta'};
all_sec = {'sec_neg10','sec_neg5','sec_0','sec_5','sec_10'};
all_contig = {'not_soz','soz'};

%% Locations
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'soz_overlap/'];
if exist('outFolder','dir') == 0
    mkdir(outFolder)
end

%% Initialize cross time comparisons
t_text = cell(length(metrics),length(all_sec),...
    length(all_freq));

% Loop over frequencies and times
for freq_idx = 1:length(all_freq)
freq = all_freq{freq_idx};

for sec_idx = 1:length(all_sec)

    
sec_text = all_sec{sec_idx};

% Initialize arrays
all_p = [];
all_z_pts = {};
all_t = [];
all_df = [];

% Initialize arrays for transitivity
all_p_trans = [];
all_t_trans = [];

% Initialize arrays for hub stability
all_p_hub = [];
all_t_hub = [];

% loop over metrics
for metric = 1:n_metrics
    
    if metric == 8
        soz_test(metric).trans.values = [];
    end
    
    % Initialize some variables
    soz_test(metric).z = [];
    soz_test(metric).rm_soz = [];
    soz_test(metric).rm_not_soz = {};
    soz_test(metric).hub.means = [];

    % loop over patients
    for i = 1:length(soz_overlap)
        if isfield(soz_overlap(i),freq) == 0, continue; end
        if isempty(soz_overlap(i).(freq)) == 1, continue; end
        
        
        
        for contig_idx = [1,2] % 1 is spare soz, 2 is remove soz
            contig_text = all_contig{contig_idx};
            if isfield(soz_overlap(i).(freq),contig_text) == 0, continue; end

            if isfield(soz_overlap(i).(freq).(contig_text),sec_text) == 0, continue; end

            % Get base
            base = soz_overlap(i).(freq).(contig_text).(sec_text);

            % Get agreement metric
            measure = base.(metrics{metric})';
            
            % Test to probe the direction of the semi-significant findings
            % for transitivity. Look at signed relative difference
            if metric == 8 % transitivity
                if contig_idx == 2
                    signed_trans_target = measure;
                elseif contig_idx == 1
                    signed_trans_spare = measure;
                end
            end

            % If it's a global metric, take absolute value and make negative
            % This is because I am just interested in the absolute
            % difference from the true value; negative to make it go in the
            % same direction as nodal measure (higher = more agreement)
            if global_metric(metric) == 1
                measure = -abs(measure);
            end
            
            % Get if hub stays same
            if metric <= 5
                same_hub = base.same_hub.(hub_metrics{metric});
                if isempty(pt(i).soz.nums) == 1
                    same_hub = nan;
                end
                if contig_idx == 2
                    soz_rm_hub_same = same_hub;
                elseif contig_idx == 1
                    not_soz_rm_hub_same = mean(same_hub);
                end
            end
            
            
            % No measure if soz empty
            if isempty(pt(i).soz.nums) == 1
                measure = nan;
                if metric == 8
                    signed_trans_target = nan;
                    signed_trans_spare = nan;
                end
            end
            
            if contig_idx == 2
                soz_rm_measure = measure;
            elseif contig_idx == 1
                not_soz_rm_measure = measure;
            end
      
        
        end
        
        if soz_rm_measure== 1
            error('what');
        end
        
        % Plot the distribution of agreement measures for the non-soz
        % removal and the agreement measure for the soz removal
        if 0 && metric == 8
            figure 
            plot(signed_trans_spare,'o')
            hold on
            xl = get(gca,'xlim');
            plot([xl(1) xl(2)],[signed_trans_target signed_trans_target])
            pause
            close(gcf)
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
        %% Take means
        % I think this is the only way to directly compare
        if exist('soz_rm_measure','var') == 0, continue; end
        if exist('not_soz_rm_measure','var') == 0, continue; end
        z = [mean(not_soz_rm_measure),soz_rm_measure];
            
        end
        
        soz_test(metric).name = metrics{metric};
        soz_test(metric).z = [soz_test(metric).z;z];
        soz_test(metric).rm_soz = [soz_test(metric).rm_soz;soz_rm_measure];
        soz_test(metric).rm_not_soz{i} = not_soz_rm_measure;
        if metric <= 5
        soz_test(metric).hub.means = [soz_test(metric).hub.means;not_soz_rm_hub_same,soz_rm_hub_same];
        end
        
        % add the transitivity signed relative differences
        if metric == 8
            soz_test(metric).trans.values = [soz_test(metric).trans.values;...
                mean(signed_trans_spare),signed_trans_target];
        end

    end
    
    if isempty(soz_test(metric).z) == 1, continue; end
    
    
    %% Paired T test on the "z scores" (in this case, the agreements)
    
    % the first column is when we spare the soz, the second column is when
    % we target the soz, and so a positive t-statistic means that the
    % agreement is higher when we spare the SOZ (so network more affected
    % when we target the soz)
    [~,p,~,stats] = ttest(soz_test(metric).z(:,1),soz_test(metric).z(:,2));
    soz_test(metric).stats.p = p;
    soz_test(metric).stats.t = stats.tstat;
    soz_test(metric).stats.df = stats.df;
    all_df = [all_df;stats.df];
    all_t = [all_t;stats.tstat];
    all_p = [all_p;p];
    all_z_pts{metric} = soz_test(metric).z(:,1)-soz_test(metric).z(:,2);
    
    %% Paired t-test on hub stability
    if metric <= 5
    
    [~,p,~,stats] = ttest(soz_test(metric).hub.means(:,1),soz_test(metric).hub.means(:,2));
    soz_test(metric).hub.stats.p = p;
    soz_test(metric).hub.stats.t = stats.tstat;
    soz_test(metric).hub.stats.df = stats.df;
    all_p_hub = [all_p_hub;p];
    all_t_hub = [all_t_hub;stats.tstat];
    end
    
    %% Paired T test for transitivity
    if metric == 8
        % The first column is spare soz and second is remove soz, and so a
        % positive t statistic means that the change in the transitivity is 
        % more positive (it goes up more) when we spare the soz. I end up
        % seeing negative t stats, meaning transitivity goes up more when
        % we target the soz, suggesting that the SOZ was acting to LOWER
        % transitivity.
        [~,p,~,stats] = ttest(soz_test(metric).trans.values(:,1),soz_test(metric).trans.values(:,2)); 
        soz_test(metric).trans.stats.p = p;
        soz_test(metric).trans.stats.t = stats.tstat;
        soz_test(metric).trans.stats.df = stats.df;
        soz_test(metric).trans.stats.means = [nanmean(soz_test(metric).trans.values(:,1)),...
            nanmean(soz_test(metric).trans.values(:,2))];
        all_p_trans = [all_p_trans;p];
        all_t_trans = [all_t_trans;stats.tstat];
        
        fprintf(['For time %s and %s frequency, mean relative difference in transitivity is\n'...
            '%1.2f when we spare the SOZ and %1.2f when we target the SOZ\n'...
            'tstat = %1.2f, p = %1.3f\n'],...
            all_sec{sec_idx},all_freq{freq_idx},nanmean(soz_test(metric).trans.values(:,1)),...
            nanmean(soz_test(metric).trans.values(:,2)),stats.tstat,p);
    end
end


%% Table
if 0
    fprintf(['For time %s and %s frequency, table of p-values and t-statistics for\n'...
        'whether agreement is higher when removing non-soz electrodes is:\n'],all_sec{sec_idx},...
        all_freq{freq_idx});

    table(metrics',all_t,all_df,all_p,'VariableNames',{'Metric','Tstat','df','p'})
end

for i = 1:size(all_t)
    t_text{i,sec_idx,freq_idx} = pretty_tstat(all_t(i),all_p(i),length(metrics));
end


for i = 1:size(all_t_trans)
    trans_text{i,sec_idx,freq_idx} = pretty_tstat(all_t_trans(i),all_p_trans(i),length(metrics));
end

for i = 1:size(all_t_hub)
    hub_text{i,sec_idx,freq_idx} = pretty_tstat(all_t_hub(i),all_p_hub(i),5);
end
%{
for i = 1:size(all_t)
    if all_p(i) < 0.001/length(metrics)
        extra = '***';
    elseif all_p(i) < 0.01/length(metrics)
        extra = '**';
    elseif all_p(i) < 0.05/length(metrics)
        extra = '*';
    else
        extra = '';
    end
    t_text{i,sec_idx,freq_idx} = sprintf('%1.2f%s',all_t(i),extra);
end
%}

%table(char(t_text(:,sec_idx,freq_idx)))


if do_plot == 1
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
    plot([i-0.3,i + 0.3],[nanmean(all_z_pts{i}),nanmean(all_z_pts{i})],'color',cols(i,:),'linewidth',3);

    xticks(1:length(all_z_pts))
    xticklabels((metric_names))
    xlim([0.7 length(all_z_pts) + 0.3])
    title({'Difference in metric agreement between seizure onset zone-sparing','and seizure onset zone-targeted resampling'})
    ylabel('Difference in agreement');
    set(gca,'fontsize',20)
    fix_xticklabels(gca,0.1,{'FontSize',20});
end
plot(get(gca,'xlim'),[0 0],'k--','linewidth',2);
for i = 1:length(stars)
    text(i + 0.15, 0.82,stars{i},'fontsize',50);
end
if save_plot == 1
    print(gcf,[outFolder,'soz_overlap_',freq,'_',sec_text],'-depsc');
end
end


end

end

fprintf('Main table:\n');
table(char(t_text(:,1,1)),char(t_text(:,2,1)),char(t_text(:,3,1)),...
    char(t_text(:,4,1)),char(t_text(:,5,1)),char(t_text(:,3,2)),'VariableNames',...
    {all_sec{1},all_sec{2},all_sec{3},all_sec{4},all_sec{5},all_freq{2}})



fprintf('Trans table:\n');
% Positive t statistics mean that the transitivity is HIGHER when we spare
% the SOZ than when we remove the SOZ
table((trans_text(:,1,1)),(trans_text(:,2,1)),(trans_text(:,3,1)),...
    (trans_text(:,4,1)),(trans_text(:,5,1)),(trans_text(:,3,2)),'VariableNames',...
    {all_sec{1},all_sec{2},all_sec{3},all_sec{4},all_sec{5},all_freq{2}})

fprintf('Hub table:\n');
% Positive t statistics mean that the hub stability is HIGHER when we spare
% the SOZ
table((hub_text(:,1,1)),(hub_text(:,2,1)),(hub_text(:,3,1)),...
    (hub_text(:,4,1)),(hub_text(:,5,1)),(hub_text(:,3,2)),'VariableNames',...
    {all_sec{1},all_sec{2},all_sec{3},all_sec{4},all_sec{5},all_freq{2}})

end

