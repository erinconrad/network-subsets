function network_stats(whichPts)

%{

This function takes an adjacency matrix A, calculates global and nodal
network measures, and then randomly resamples the network, retaining
specified fractions of electrodes, and then gets statistics on how the
network measures change.


Stuff to add:
- get other network metrics 
    - path length
    (https://www.sciencedirect.com/science/article/pii/S1388245715012584)
    - clustering coefficient
    (https://www.sciencedirect.com/science/article/pii/S1388245715012584)
    - node heterogeneity
    - epileptor model - The Virtual Brain

- take number N of resected electrodes and randomly move them around so
still N contiguous electrodes and recalculate the control centrality of
the "resected region"

- which adjacency matrices to use?
    - 5 seconds before
    - right at seizure onset
    - 10 seconds after
- aim3/results/patientID/aim3/multiband - minutes long, nchxnchxTxfband
- multiple freq bands
    - start with high gamma

- brain connectivity toolbox


%}

tic

%% Parameters

% 1 if doing SOZ analysis, 0 if doing main analysis
do_soz_analysis = 0;

doPlots = 1;

% add to existing stats array? Or re-write with new stats array?
merge = 1;

% How many random resamples to do of each fraction
n_perm = 1e2;

% What fraction of nodes to retain
if do_soz_analysis == 1
    e_f = 0.8;
else
    e_f = [0.2 0.4 0.6 0.8 1];
end
n_f = length(e_f);

% Which freq?
freq = 'high_gamma';


%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);

%% Load the output structure to add more info to it
if merge == 1
    if do_soz_analysis == 1
        if exist([resultsFolder,'basic_metrics/soz.mat'],'file') ~= 0
            load([resultsFolder,'basic_metrics/soz.mat']);
        else
            soz = struct;
        end
    else
        if exist([resultsFolder,'basic_metrics/stats.mat'],'file') ~= 0
            load([resultsFolder,'basic_metrics/stats.mat']);
        else
            stats = struct;
        end
    end
else
    if do_soz_analysis == 1
        soz = struct;
    else
        stats = struct;
    end
end


%% Loop through patients, times, and whether contig or random electrodes
for which_sec = [-5 0] % 0 means start time of the seizure, -5 is 5 seconds before
for contig = [1 0] %1 means semi-contiguous set of electrodes, 0 means random electrodes

% Loop through patients
for whichPt = whichPts

    % Make result folder
    name = pt(whichPt).name;
    
    fprintf('Doing %s\n',name);
    
    outFolder = [resultsFolder,'basic_metrics/',name,'/'];
    if exist(outFolder,'dir') == 0
        mkdir(outFolder);
    end

    if contig == 1
        contig_text = 'contiguous';
    else
        contig_text = 'random';
    end
    
    if which_sec < 0
        sec_text = sprintf('sec_neg%d',abs(which_sec));
    else
        sec_text = sprintf('sec_%d',abs(which_sec));
    end
    
    % Continue if we've already done it
    if merge == 1
        if do_soz_analysis == 1
            stats = soz;
        end
        if length(stats) >= whichPt
            if isfield(stats(whichPt),'eff') == 1
                if isfield(stats(whichPt).eff,(contig_text)) == 1
                    if isfield(stats(whichPt).eff.(contig_text),(sec_text)) == 1
                        if isfield(stats(whichPt).eff.(contig_text).(sec_text),'true') == 1
                            if isempty(stats(whichPt).eff.(contig_text).(sec_text).true) == 0
                                fprintf('Did %s, skipping\n',name);
                                continue
                            end
                        end
                    end
                end
            end 
        end
    end

    % Get adjacency matrix
    [adj,~] = reconcileAdj(pt,whichPt);

    if strcmp(freq,'high_gamma') == 1
        A_all = adj(4).data;
        if contains(adj(4).name,'highgamma') == 0
            error('This isn''t gamma!'\n');
        end
    end

    % Start with the middle and add which second
    A = squeeze(A_all(ceil(size(A_all,1)/2)+which_sec,:,:));

    %% Get true control centrality
    c_c = control_centrality(A);
    fprintf('There are %d synchronizing and %d desynchronizing nodes.\n',...
        sum(c_c<0),sum(c_c>0));
    
    % Get identity of node with lowest control centrality
    [~,min_cc_true] = min(c_c);
    
    % Get locations
    locs = pt(whichPt).new_elecs.locs; % all electrode locations
    
    %% Get regional control centrality
    
    if isempty(pt(whichPt).resec) == 0
        
        % Get number of resected electrodes
        num_resec = length(pt(whichPt).resec.nums);

        % Get regional control centralities
        [cc_regional,elecs_regional] = regional_control_centrality(A,num_resec,locs,1);

        % Get identity of region with lowest regional control centrality
        [~,min_cc_regional_true] = min(cc_regional);

        % Get electrodes in region with lowest cc
        elecs_regional_min = elecs_regional(min_cc_regional_true,:);
    else
        % Not doing it if we don't have resection data
        elecs_regional_min = nan;
    end

    %% Get true synchronizability
    sync = synchronizability(A);

    %% Get true betweenness centrality
    bc = betweenness_centrality(A,1);

    %% Get true node strength
    ns = node_strength(A);

    %% Get true global efficiency
    eff = efficiency_wei(A, 0);

    %% Resample network and get new metrics
    % all_c_c is nch x n_f x n_perm size matrix
    [all_c_c,all_ns,all_bc,all_sync,all_eff,overlap_soz,dist_soz,...
        overlap_resec,dist_resec,temp_centroid_min] = ...
        resampleNetwork(A,n_perm,e_f,contig,pt,whichPt,adj);

    %% Initialize SMC and rho arrays for node-level metrics

    % Control centrality stuff
    SMC_cc = zeros(n_f,n_perm);
    rho_cc = zeros(n_f,n_perm);
    true_cc_most_sync = zeros(n_f,n_perm);
    min_cc_resample_loc = zeros(n_f,n_perm,3);

    % Betweenness centrality stuff
    rho_bc = zeros(n_f,n_perm);

    % Node strength centrality stuff
    rho_ns = zeros(n_f,n_perm);

    %% Loop over each fraction and get various stats
    for f = 1:n_f
        c_c_f = squeeze(all_c_c(:,f,:));
        ns_f = squeeze(all_ns(:,f,:));
        bc_f = squeeze(all_bc(:,f,:));
        
        % Loop over each permutation
        for i_p = 1:n_perm
            c_c_f_p = squeeze(c_c_f(:,i_p));
            ns_f_p = squeeze(ns_f(:,i_p));
            bc_f_p = squeeze(bc_f(:,i_p));

            %% Do stats on control centrality

            % Get rho and SMC
            [rho_cc(f,i_p),SMC_cc(f,i_p)] = doStats(c_c,c_c_f_p);

            if sum(isnan(c_c_f_p)) == length(c_c_f_p)
                true_cc_most_sync(f,i_p) = nan;
            else

                % Get the identity of the most synchronizing node (the one we would
                % tell the surgeons to resect)
                [~,ch_most_sync] = min(c_c_f_p);
                
                % Get distances between lowest control centrality electrode
                % in resampled network and original network
                min_cc_resample_loc(f,i_p,:) = locs(ch_most_sync,:);

                % Fill up SMC and rho arrays
                true_cc_most_sync(f,i_p) = c_c(ch_most_sync);
                
            end


            %% Do Spearman rank on node strength and betweenness centrality
            % For these, they are always non-negative and so SMC doesn't make
            % sense
            [rho_ns(f,i_p),~] = doStats(ns,ns_f_p);
            [rho_bc(f,i_p),~] = doStats(bc,bc_f_p);


        end

    end

    %% Agreement measures, which I will include in both analyses
    % Global measures of agreement: relative difference
    rel_sync = (all_sync-sync)/sync;
    rel_eff = (all_eff-eff)/eff;

    % Nodal measures: "average" the SRCs. To average the SRCs, I
    % apply a Fisher's transformation, average the z values, and then back
    % transform to get an average rho.
    rho_mean_cc = average_rho(rho_cc,2);
    rho_mean_ns = average_rho(rho_ns,2);
    rho_mean_bc = average_rho(rho_bc,2);
     
    %% Fill up stats structures
    
    if do_soz_analysis == 0
    
        %% Doing the main analysis (not dependence on distance from SOZ/resection zone)
        % For this, we need measures of variability AND agreement
       
        %% Nodal measure of variability
        % Nodal measure: relative std (std across permuations divided by
        % std across electrodes)
        cc_rel_std = rel_std_nodal(all_c_c,c_c);
        ns_rel_std = rel_std_nodal(all_ns,ns);
        bc_rel_std = rel_std_nodal(all_bc,bc);
        
        % Global measures: we will do relative std (std across
        % permutations divided by std across patients). And so we will
        % record the standard deviation for now.
        
        % How often would we resect the wrong piece of brain?
        resect_wrong = sum((true_cc_most_sync > 0),2)/n_perm;
        
        % Mean and std of location of lowest cc in resampled network
        cc_res_mean = squeeze(nanmean(min_cc_resample_loc,2));
        cc_res_std = squeeze(nanstd(min_cc_resample_loc,0,2));

        % Mean and std of location of centroid of lowest cc region in resampled network
        cc_region_res_mean = squeeze(nanmean(temp_centroid_min,2));
        cc_region_res_std = squeeze(nanstd(temp_centroid_min,0,2));
        
        stats(whichPt).name = name;

        % control centrality
        stats(whichPt).(contig_text).(sec_text).cc.true = c_c;
        stats(whichPt).(contig_text).(sec_text).cc.rel_std = cc_rel_std;
        stats(whichPt).(contig_text).(sec_text).cc.resect_wrong = resect_wrong;
        stats(whichPt).(contig_text).(sec_text).cc.rho_mean = rho_mean_cc;
        
        % Mean and STD of resampled min cc
        stats(whichPt).(contig_text).(sec_text).min_cc.res_mean = cc_res_mean;
        stats(whichPt).(contig_text).(sec_text).min_cc.res_std = cc_res_std;
        stats(whichPt).(contig_text).(sec_text).min_cc.true = min_cc_true;
        
        % Regional cc
        stats(whichPt).(contig_text).(sec_text).regional_cc.true = elecs_regional_min;
        stats(whichPt).(contig_text).(sec_text).regional_cc.res_mean = cc_region_res_mean;
        stats(whichPt).(contig_text).(sec_text).regional_cc.res_std = cc_region_res_std;

        % node strength
        stats(whichPt).(contig_text).(sec_text).ns.rel_std = ns_rel_std;
        stats(whichPt).(contig_text).(sec_text).ns.rho_mean = rho_mean_ns;

        % betweenness centrality
        stats(whichPt).(contig_text).(sec_text).bc.rel_std = bc_rel_std;
        stats(whichPt).(contig_text).(sec_text).bc.rho_mean = rho_mean_bc;

        % synchronizability
        stats(whichPt).(contig_text).(sec_text).sync.std = std(all_sync,0,2);
        stats(whichPt).(contig_text).(sec_text).sync.true = sync;
        stats(whichPt).(contig_text).(sec_text).sync.rel_diff = rel_sync;
        

        % efficiency
        stats(whichPt).(contig_text).(sec_text).eff.std = std(all_eff,0,2);
        stats(whichPt).(contig_text).(sec_text).eff.true = eff;
        stats(whichPt).(contig_text).(sec_text).eff.rel_diff = rel_eff;

        save([resultsFolder,'basic_metrics/stats.mat'],'stats');
    elseif do_soz_analysis == 1
        %% Do the analysis of dependence of agreement on distance from important things
        
        
        
        soz(whichPt).(contig_text).(sec_text).rho_cc = rho_mean_cc;
        soz(whichPt).(contig_text).(sec_text).rho_bc = rho_mean_bc;
        soz(whichPt).(contig_text).(sec_text).rho_ns = rho_mean_ns;
        soz(whichPt).(contig_text).(sec_text).sync = rel_sync;
        soz(whichPt).(contig_text).(sec_text).eff = rel_eff;
        soz(whichPt).(contig_text).(sec_text).dist_soz = dist_soz;
        soz(whichPt).(contig_text).(sec_text).dist_resec = dist_resec;
        soz(whichPt).(contig_text).(sec_text).overlap_soz = overlap_soz;
        soz(whichPt).(contig_text).(sec_text).overlap_resec = overlap_resec;
        
        save([resultsFolder,'basic_metrics/soz.mat'],'soz');
    end

    %% Plots
  
    if doPlots == 1
        
        if do_soz_analysis == 1
        
        for i = 1:n_f
            %% Plot metrics as a function of relation to SOZ
            figure
            set(gcf,'Position',[92 126 1349 672]);
            % Distance from SOZ, global metrics
            subplot(2,2,1)
            scatter(dist_soz(i,:),rel_sync(i,:),100,'filled','b');
            hold on
            scatter(dist_soz(i,:),rel_eff(i,:),100,'filled','r');
         %   legend('Synchronizability','True synchronizability','Global efficiency','True global efficiency',...
         %       'location','northeastoutside');
            title('Global metrics as a function of distance from SOZ');
            xlabel('Distance from SOZ');
            ylabel('Metric');
            set(gca,'fontsize',20);

            % Distance from SOZ, nodal metrics
            subplot(2,2,2)
            scatter(dist_soz(i,:),rho_cc(i,:),100,'filled','b');
            hold on
            scatter(dist_soz(i,:),rho_bc(i,:),100,'filled','r');
            scatter(dist_soz(i,:),rho_ns(i,:),100,'filled','g');
            %legend('Control centrality','Betweenness centrality','Node strength','location','northeastoutside');
            title('Correlation of nodal metrics');
            xlabel('Distance from SOZ');
            ylabel('Correlation of metric');
            set(gca,'fontsize',20);

            % Overlap with SOZ, global metrics
            subplot(2,2,3)
            scatter(overlap_soz(i,:)*100,rel_sync(i,:),100,'filled','b');
            hold on
            scatter(overlap_soz(i,:)*100,rel_eff(i,:),100,'filled','r');
            %legend('Synchronizability','True synchronizability','Global efficiency','True global efficiency');
            title('Global metrics as a function of overlap with SOZ');
            xlabel('% of electrodes removed that were in SOZ');
            ylabel('Metric');
            set(gca,'fontsize',20);

            % Overlap with SOZ, nodal metrics
            subplot(2,2,4)
            scatter(overlap_soz(i,:)*100,rho_cc(i,:),100,'filled','b');
            hold on
            scatter(overlap_soz(i,:)*100,rho_bc(i,:),100,'filled','r');
            scatter(overlap_soz(i,:)*100,rho_ns(i,:),100,'filled','g');
            %legend('Control centrality','Betweenness centrality','Node strength');
            title('Correlation of nodal metrics');
            xlabel('% of electrodes removed that were in SOZ');
            ylabel('Correlation of metric');
            set(gca,'fontsize',20);
            print(gcf,[outFolder,'soz_',contig_text,sec_text],'-depsc');
            close(gcf)

            %% Plot metrics as a function of relation to resection
            figure
            set(gcf,'Position',[92 126 1349 672]);
            % Distance from resec, global metrics
            subplot(2,2,1)
            scatter(dist_resec(i,:),rel_sync(i,:),100,'filled','b');
            hold on
            scatter(dist_resec(i,:),rel_eff(i,:),100,'filled','r');
           % legend('Synchronizability','True synchronizability','Global efficiency','True global efficiency',...
           %     'location','northeastoutside');
            title('Global metrics');
            xlabel('Distance from resected electrodes');
            ylabel('Metric');
            set(gca,'fontsize',20);

            % Distance from resec, nodal metrics
            subplot(2,2,2)
            scatter(dist_resec(i,:),rho_cc(i,:),100,'filled','b');
            hold on
            scatter(dist_resec(i,:),rho_bc(i,:),100,'filled','r');
            scatter(dist_resec(i,:),rho_ns(i,:),100,'filled','g');
           % legend('Control centrality','Betweenness centrality','Node strength',...
           %     'location','northeastoutside');
            title('Correlation of nodal metrics');
            xlabel('Distance from resected electrodes');
            ylabel('Correlation of metric');
            set(gca,'fontsize',20);

            % Overlap with resec, global metrics
            subplot(2,2,3)
            scatter(overlap_resec(i,:)*100,rel_sync(i,:),100,'filled','b');
            hold on
            scatter(overlap_resec(i,:)*100,rel_eff(i,:),100,'filled','r');
           % legend('Synchronizability','True synchronizability','Global efficiency','True global efficiency');
            title('Global metrics');
            xlabel('% of electrodes removed that were in resected region');
            ylabel('Metric');
            set(gca,'fontsize',20);

            % Overlap with resec, nodal metrics
            subplot(2,2,4)
            scatter(overlap_resec(i,:)*100,rho_cc(i,:),100,'filled','b');
            hold on
            scatter(overlap_resec(i,:)*100,rho_bc(i,:),100,'filled','r');
            scatter(overlap_resec(i,:)*100,rho_ns(i,:),100,'filled','g');
           % legend('Control centrality','Betweenness centrality','Node strength');
            title('Correlation of nodal metrics');
            xlabel('% of electrodes removed that were in resected region');
            ylabel('Correlation');
            set(gca,'fontsize',20);
            print(gcf,[outFolder,'resec_',contig_text,sec_text],'-depsc');
            close(gcf)
        
        end
        else
        
        %% Plot variability metrics as a function network retained
        figure
        set(gcf,'Position',[50 133 1211 665]);
        subplot(1,2,1)
        scatter(e_f,cc_rel_std,100,'filled');
        hold on
        scatter(e_f,bc_rel_std,100,'filled');
        scatter(e_f,ns_rel_std,100,'filled');
        ylabel('Relative std of metric')
        xlabel('Fraction of nodes retained')
        title('Metric variability by subsample size');
        legend('Control centrality','Betweenness centrality','Node strength')
        set(gca,'fontsize',20);
        
        subplot(1,2,2)
        scatter(e_f,rho_mean_cc,100,'filled');
        hold on
        scatter(e_f,rho_mean_bc,100,'filled');
        scatter(e_f,rho_mean_ns,100,'filled');
        ylabel('Spearman rank correlation of metric with original')
        xlabel('Fraction of nodes retained')
        title('Average metric agreement by subsample size');
        legend('Control centrality','Betweenness centrality','Node strength',...
            'Location','southeast')
        set(gca,'fontsize',20);
        
        
        print(gcf,[outFolder,'nodal_metrics_',contig_text,sec_text],'-depsc');
        close(gcf)
     
        end
    end
    
end

end
end


toc

end

function rel_std = rel_std_nodal(perm_metric,true_metric)
    
    % The average across electrodes of the std across permutations
    std_across_perm = nanmean(squeeze(nanstd(perm_metric,0,3)),1);
        
    % Std across electrodes of the true cc
    std_across_ch = nanstd(true_metric,0,1);

    % Relative std
    rel_std = std_across_perm./std_across_ch;

end