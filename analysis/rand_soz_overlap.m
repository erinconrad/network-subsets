function rand_soz_overlap(soz_overlap,pt)

%{
This script compares network statistics when we resample by targeting the
SOZ as opposed to removing a random group of electrodes equal in number to
the SOZ. (unlike soz_overlap_analysis, which specifically compares
targeting vs sparing the soz).
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
all_contig = {'random','soz'};

%% Locations
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'soz_overlap_random/'];
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


% loop over metrics
for metric = 1:n_metrics
    
    
    % Initialize some variables
    soz_test(metric).z = [];
    soz_test(metric).p = [];

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
            
            
            % No measure if soz empty
            if isempty(pt(i).soz.nums) == 1
                measure = nan;
            end
            
            if contig_idx == 2
                soz_rm_measure = measure;
            elseif contig_idx == 1
                random_rm_measure = measure;
            end
        
        end
        
        if soz_rm_measure== 1
            error('what');
        end
        
        % Plot the distribution of agreement measures for the random
        % removal and the agreement measure for the soz removal
        if 0 && metric == 8
            figure 
            plot(random_rm_measure,'o')
            hold on
            xl = get(gca,'xlim');
            plot([xl(1) xl(2)],[soz_rm_measure soz_rm_measure])
            pause
            close(gcf)
        end
        
        %% Compare the soz targeted removal to random removal
        
        % Compute percentile/100 - 0.5 (for a signed test statistic)
        pcntile = sum(random_rm_measure<soz_rm_measure)/length(random_rm_measure) - 0.5;
        
        % do plot
        if 0
            figure
            histogram(random_rm_measure)
            hold on
            yl = get(gca,'ylim');
            plot([soz_rm_measure soz_rm_measure],[yl(1) yl(2)])
            title(sprintf('prcntile = %1.1f, p = %1.3f',pcntile,p))
            pause
            close(gcf)
        end
        
        
        %% Fill structure
        if isempty(pt(i).soz.nums) == 1
            pcntile = nan;
        end
        soz_test(metric).name = metrics{metric};
        soz_test(metric).z = [soz_test(metric).z;pcntile];

    end
    
    if isempty(soz_test(metric).z) == 1, continue; end
    
    
    %% One sample T test on the "z scores" (in this case, the scaled and shifted percentiles)
    [~,p,~,stats] = ttest(soz_test(metric).z(~isnan(soz_test(metric).z)));
    soz_test(metric).stats.p = p;
    soz_test(metric).stats.t = stats.tstat;
    soz_test(metric).stats.df = stats.df;
    soz_test(metric).stats.all_percentiles = (soz_test(metric).z);
    
    % Add to array for plotting
    all_df = [all_df;stats.df];
    all_t = [all_t;stats.tstat];
    all_p = [all_p;p];
    all_z_pts{metric} = (soz_test(metric).z+0.5)*100;
    
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
    title({'Percentile agreement when seizure onset zone','is targeted for resampling'})
    ylabel('Percentile agreement');
    set(gca,'fontsize',20)
    fix_xticklabels(gca,0.1,{'FontSize',20});
end
plot(get(gca,'xlim'),[50 50],'k--','linewidth',2);
for i = 1:length(stars)
    text(i + 0.15, 90,stars{i},'fontsize',50);
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


end

