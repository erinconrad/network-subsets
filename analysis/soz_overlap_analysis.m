function soz_overlap_analysis(soz_overlap,pt,example,only_ilae1,do_resec)


%{
This is the analysis we added for Network Neuroscience. It compares metrics
when we target vs spare the SOZ for removal
%}

%% Parameters
rm_non_independent_pts = 0;
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

if example == 1
    all_freq = {'high_gamma'};
    all_sec = {'sec_0'};
else
    all_freq = {'high_gamma','beta'};
    all_sec = {'sec_neg10','sec_neg5','sec_0','sec_5','sec_10'};
end
if do_resec == 0
    all_contig = {'not_soz','soz'};
elseif do_resec == 1
    all_contig = {'not_resec','resec'};
end

%% Locations
if example ~= 1
    [electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
    pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
    outFolder = [resultsFolder,'soz_overlap/'];
    if exist('outFolder','dir') == 0
        mkdir(outFolder)
    end
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

all_p_dist = [];
all_z_pts = {};
all_t_dist = [];
all_df_dist = [];

% Initialize arrays for transitivity
all_p_trans = [];
all_t_trans = [];
all_p_trans_dist = [];
all_t_trans_dist = [];

% Initialize arrays for hub stability
all_p_hub = [];
all_t_hub = [];

% Initialize agreement array for example patient
agreement_ex = nan(8,2);

% loop over metrics
for metric = 1:n_metrics
    
    if metric == 8
        soz_test(metric).trans.values = [];
        soz_test(metric).trans.all_spare = [];
    end
    
    % Initialize some variables
    soz_test(metric).z = [];
    soz_test(metric).rm_soz = [];
    soz_test(metric).rm_not_soz = {};
    soz_test(metric).hub.means = [];
    soz_test(metric).spare_lower_per = [];

    % loop over patients
    for i = 1:length(soz_overlap)
        if isfield(soz_overlap(i),freq) == 0, continue; end
        if isempty(soz_overlap(i).(freq)) == 1, continue; end
        
        % Remove non ilae 1 patients
            if only_ilae1 == 1
                % get outcome
                outcome = get_ilae(pt(i).name);
                
                if outcome ~=1, continue; end
            end
        
        % Remove non-independent patients
        if rm_non_independent_pts == 1
            if ismember(i,[1 8 21 23 25 29]) == 1, continue; end
            
        end
        
        
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
        
        
        %% Get distribution of soz_removing
        % Comparing single rho from when we only remove soz to all rhos
        % when we randomly remove something that is not the soz

        % The reason this is problematic is because the group of
        % not_soz_rm_measures are not independent. For instance, there is
        % one patient where the soz electrodes comprise more than half of
        % the electrodes. And so I just take the remainder of the
        % electrodes to remove for not_soz_rm_measure. And I take that same
        % set 1,000 times. 

        if exist('soz_rm_measure','var') == 0, continue; end

        % Get the percentage of sparing agreements lower than the
        % targeted agreement
        spare_lower_per = sum(not_soz_rm_measure < soz_rm_measure)/...
            length(not_soz_rm_measure);


        %[p,h,stats] = ranksum(soz_rm_measure,not_soz_rm_measure);

        % Get the z-score
        %z = stats.zval;


        %% Take means
        % I think this is the only way to directly compare
        if exist('soz_rm_measure','var') == 0, continue; end
        if exist('not_soz_rm_measure','var') == 0, continue; end
        z = [mean(not_soz_rm_measure),soz_rm_measure];

        if example == 1
            agreement_ex(metric,:) = z;
            continue
        end



        soz_test(metric).name = metrics{metric};
        soz_test(metric).z = [soz_test(metric).z;z];
        soz_test(metric).rm_soz = [soz_test(metric).rm_soz;soz_rm_measure];
        soz_test(metric).rm_not_soz{i} = not_soz_rm_measure;
        soz_test(metric).spare_lower_per = [soz_test(metric).spare_lower_per;spare_lower_per];
        if metric <= 5
        soz_test(metric).hub.means = [soz_test(metric).hub.means;not_soz_rm_hub_same,soz_rm_hub_same];
        end

        % add the transitivity signed relative differences
        if metric == 8
            soz_test(metric).trans.values = [soz_test(metric).trans.values;...
                mean(signed_trans_spare),signed_trans_target];
            soz_test(metric).trans.all_spare = [soz_test(metric).trans.all_spare;...
                sum(signed_trans_spare<signed_trans_target)/length(signed_trans_spare)];
        end
        
        

    end
    
    if example == 1, continue; end
    
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
    
    % Do a test to see if the percentage of SOZ-sparing agreements less
    % than SOZ-targeted agreements is significantly less than 50% (I would
    % expect, by chance, that half of the patients would have >50% and half
    % would have <50%. A positive t statistic means that SOZ-sparing
    % agreements tend to be higher.
    [~,p_dist,~,stats_dist] = ttest(0.5-soz_test(metric).spare_lower_per);
    %[p_dist,~,stats_dist] = signrank(0.5-soz_test(metric).spare_lower_per);
    soz_test(metric).stats_dist.p = p_dist;
    soz_test(metric).stats_dist.t = stats_dist.tstat;
    %soz_test(metric).stats_dist.t = stats_dist.signedrank;
    soz_test(metric).stats_dist.df = stats_dist.df;
    %soz_test(metric).stats_dist.df = nan;
    all_df_dist = [all_df_dist;soz_test(metric).stats_dist.df];
    all_t_dist = [all_t_dist;soz_test(metric).stats_dist.t];
    all_p_dist = [all_p_dist;p_dist];
    
    %if metric== 8, error('look\n'); end
    
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
        
        % Do test on distribution
        [~,p_dist,~,stats_dist] = ttest(0.5-soz_test(metric).trans.all_spare);
        soz_test(metric).trans.stats_dist.p = p_dist;
        soz_test(metric).trans.stats_dist.t = stats_dist.tstat;
        soz_test(metric).trans.stats_dist.df = stats_dist.df;
        all_p_trans_dist = [all_p_trans_dist;p_dist];
        all_t_trans_dist = [all_t_trans_dist;stats_dist.tstat];
        
    end
end

if example == 1
    % plot figure for example patient
    figure
    set(gcf,'position',[1         416        1440         382]);
    ha(1) = subplot(1,2,1);
    set(ha(1),'position',[0.05 0.07 0.48 0.85])
    plot(0.85:2:8.85,agreement_ex(1:5,1),'o','markersize',15,'MarkerFaceColor','auto');
    hold on
    plot(1.15:2:9.15,agreement_ex(1:5,2),'o','markersize',15,'MarkerFaceColor','auto');
    xticks(1:2:9)
    xticklabels(metric_names(1:5))
    ylabel('Resampled-original metric agreement');
    legend({'SOZ-sparing','SOZ-targeting'},'location','southeast')
    title('Nodal metrics')
    set(gca,'fontsize',13)
    ha(2) = subplot(1,2,2);
    set(ha(2),'position',[0.58 0.07 0.42 0.85])
    plot(0.85:1:2.85,agreement_ex(6:8,1),'o','markersize',15,'MarkerFaceColor','auto');
    hold on
    plot(1.15:1:3.15,agreement_ex(6:8,2),'o','markersize',15,'MarkerFaceColor','auto');
    legend({'SOZ-sparing','SOZ-targeting'},'location','southeast')
    xticks(1:1:3)
    xticklabels(metric_names(6:8))
    title('Global metrics')
    set(gca,'fontsize',13)
    return
end


%% Table
if 0
    fprintf(['For time %s and %s frequency, table of p-values and t-statistics for\n'...
        'whether agreement is higher when removing non-soz electrodes is:\n'],all_sec{sec_idx},...
        all_freq{freq_idx});

    table(metrics',all_t,all_df,all_p,'VariableNames',{'Metric','Tstat','df','p'})
end

% T stats
for i = 1:size(all_t)
    t_text{i,sec_idx,freq_idx} = pretty_tstat(all_t(i),all_p(i),length(metrics));
end

% T stats for distribution
for i = 1:size(all_t_dist)
    t_text_dist{i,sec_idx,freq_idx} = pretty_tstat(all_t_dist(i),all_p_dist(i),length(metrics));
end

% Signed transitivity
for i = 1:size(all_t_trans)
    trans_text{i,sec_idx,freq_idx} = pretty_tstat(all_t_trans(i),all_p_trans(i),length(metrics));
    trans_text_dist{i,sec_idx,freq_idx} = pretty_tstat(all_t_trans_dist(i),all_p_trans_dist(i),length(metrics));
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
% Main t stats
table(char(t_text(:,1,1)),char(t_text(:,2,1)),char(t_text(:,3,1)),...
    char(t_text(:,4,1)),char(t_text(:,5,1)),char(t_text(:,3,2)),'VariableNames',...
    {all_sec{1},all_sec{2},all_sec{3},all_sec{4},all_sec{5},all_freq{2}})


% T stats for distribution
fprintf('Distribution table:\n');
table(char(t_text_dist(:,1,1)),char(t_text_dist(:,2,1)),char(t_text_dist(:,3,1)),...
    char(t_text_dist(:,4,1)),char(t_text_dist(:,5,1)),char(t_text_dist(:,3,2)),'VariableNames',...
    {all_sec{1},all_sec{2},all_sec{3},all_sec{4},all_sec{5},all_freq{2}})

fprintf('Trans table:\n');
% Positive t statistics mean that the transitivity is HIGHER when we spare
% the SOZ than when we remove the SOZ (signed transitivity)
table((trans_text(:,1,1)),(trans_text(:,2,1)),(trans_text(:,3,1)),...
    (trans_text(:,4,1)),(trans_text(:,5,1)),(trans_text(:,3,2)),'VariableNames',...
    {all_sec{1},all_sec{2},all_sec{3},all_sec{4},all_sec{5},all_freq{2}})

fprintf('Distribution trans table:\n');
table((trans_text_dist(:,1,1)),(trans_text_dist(:,2,1)),(trans_text_dist(:,3,1)),...
    (trans_text_dist(:,4,1)),(trans_text_dist(:,5,1)),(trans_text_dist(:,3,2)),'VariableNames',...
    {all_sec{1},all_sec{2},all_sec{3},all_sec{4},all_sec{5},all_freq{2}})

fprintf('Hub table:\n');
% Positive t statistics mean that the hub stability is HIGHER when we spare
% the SOZ
table((hub_text(:,1,1)),(hub_text(:,2,1)),(hub_text(:,3,1)),...
    (hub_text(:,4,1)),(hub_text(:,5,1)),(hub_text(:,3,2)),'VariableNames',...
    {all_sec{1},all_sec{2},all_sec{3},all_sec{4},all_sec{5},all_freq{2}})

end

