function network_stats(whichPts,do_soz_analysis)

% do 100 random and average for fake matrix

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

% do_soz_analysis: 1 if doing SOZ analysis, 0 if doing main analysis

doSave = 0;
doPlots = 1;

% add to existing stats array? Or re-write with new stats array?
merge = 1;

% What fraction of nodes to retain
if do_soz_analysis == 1
    e_f = 0.8;
    contigs = 1;
else
    e_f = 0.8;%[0.2 0.4 0.6 0.8 1];
    contigs = [0 1];
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
for contig = contigs %1 means semi-contiguous set of electrodes, 0 means random electrodes

% Loop through patients
for whichPt = whichPts
    
    if unique(pt(whichPt).new_elecs.locs) == -1
        continue
    end
    
    if do_soz_analysis == 1
        % Here, not taking random samples, but rather systematically going
        % through all contiguous chunks
        n_perm = length(pt(whichPt).new_elecs.electrodes);
    else
        n_perm = 1e2;
    end

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
            if isfield(stats(whichPt),(contig_text)) == 1
                if isfield(stats(whichPt).(contig_text),(sec_text)) == 1
                    
                    fprintf('Did %s, skipping\n',name);
                   % continue
                end
            end 
        end
    end

    %% Get adjacency matrix
    [adj,~] = reconcileAdj(pt,whichPt);

    if strcmp(freq,'high_gamma') == 1
        A_all = adj(4).data;
        if contains(adj(4).name,'highgamma') == 0
            error('This isn''t gamma!'\n');
        end
    end

    % Start with the middle and add which second
    A = squeeze(A_all(ceil(size(A_all,1)/2)+which_sec,:,:));

    %% Get true metrics
    fprintf('Getting true metrics\n');
    
    %% Control centrality
    c_c = control_centrality(A);
    
    % Get identity of node with lowest control centrality
    [~,min_cc_true] = min(c_c);
    
    % Get locations
    locs = pt(whichPt).new_elecs.locs; % all electrode locations
    
    %% Regional control centrality
    
    if isempty(pt(whichPt).resec) == 0
        
        % Get number of resected electrodes
        num_resec = length(pt(whichPt).resec.nums);

        % Get regional control centralities
        [cc_regional,elecs_regional] = regional_control_centrality(A,num_resec,locs,1);

        % Get identity of region with lowest regional control centrality
        [~,min_cc_regional_true] = min(cc_regional);

        % Get electrodes in region with lowest cc
        elecs_regional_min = elecs_regional(min_cc_regional_true,:);
        
        % Get dice score to compare n electrodes in min cc region and
        % electrodes that themselves have the n smallest cc regions (should
        % be similar but probably not equivalent)
        %{
        [cc_regional_sorted,I] = sort(cc_regional);
        is_min_cc_1 = zeros(size(locs,1),1);
        is_min_cc_1(I(1:num_resec)) = 1;
        
        is_min_cc_2 = zeros(size(locs,1),1);
        is_min_cc_2(elecs_regional_min) = 1;
        ds = dice_score(is_min_cc_1,is_min_cc_2);
        %}
        
    else
        % Not doing it if we don't have resection data
        elecs_regional_min = nan;
        cc_regional = nan;
    end

    %% Get true synchronizability
    sync = synchronizability(A);
    sync_fake = synchronizability(generate_fake_graph(A));
    sync_norm = sync/sync_fake;

    %% Get true betweenness centrality
    bc = betweenness_centrality(A,1);
    
    %% Get true eigenvector centrality
    ec = eigenvector_centrality_und(A);
    
    %% Get true clustering coefficient
    clust = clustering_coef_wu(A);
    
    %% Get local efficiency
   % le = efficiency_wei(A,1);

    %% Get true node strength
    ns = node_strength(A);
    
    %% Get true participation coefficient
    [Ci,~] = modularity_und(A);
    par = participation_coef(A,Ci);

    %% Get true global efficiency and efficiency for fake network
    eff = efficiency_wei(A, 0);
    eff_fake = efficiency_wei(generate_fake_graph(A),0);
    eff_norm = eff/eff_fake;
    
    %% Get true transitivity
    trans = transitivity_wu(A);
    trans_norm = trans/...
            transitivity_wu(generate_fake_graph(A));

    fprintf('Got true metrics, now resampling network...\n');
    %% Resample network and get new metrics
    % all_c_c is nch x n_f x n_perm size matrix
    [all_c_c,all_ns,all_bc,all_sync,all_eff,overlap_soz,dist_soz,...
        overlap_resec,dist_resec,elecs_min,...
        all_par,all_trans,avg_par_removed,avg_bc_removed,...
        all_sync_norm,all_eff_norm,all_trans_norm,all_ec,...
        all_clust,all_le,cc_reg] = ...
        resampleNetwork(A,n_perm,e_f,contig,pt,whichPt,adj,do_soz_analysis);

    %% Initialize SMC and rho arrays for node-level metrics

    % Control centrality stuff
    rho_cc = zeros(n_f,n_perm);
    most_sync = zeros(n_f,n_perm);

    % Betweenness centrality stuff
    rho_bc = zeros(n_f,n_perm);

    % Node strength centrality stuff
    rho_ns = zeros(n_f,n_perm);
    
    % Participation coefficient stuff
    rho_par = zeros(n_f,n_perm);
    
    % Eigenvalue centality stuff
    rho_ec = zeros(n_f,n_perm);
    
    % Local efficiency stuff
    rho_le = zeros(n_f,n_perm);
    
    % Clustering coefficient stuff
    rho_clust = zeros(n_f,n_perm);
    
    % Regional control centrality
    rho_cc_reg = zeros(n_f,n_perm);

    % Rhos for things in the resection zone
    rho_cc_resec = zeros(n_f,n_perm);
    rho_bc_resec = zeros(n_f,n_perm);
    rho_ns_resec = zeros(n_f,n_perm);
    rho_par_resec = zeros(n_f,n_perm);
    rho_ec_resec = zeros(n_f,n_perm);
    rho_clust_resec = zeros(n_f,n_perm);
    rho_le_resec = zeros(n_f,n_perm);
    
    %% Loop over each fraction and get various stats
    for f = 1:n_f

        
        c_c_f = squeeze(all_c_c(:,f,:));
        ns_f = squeeze(all_ns(:,f,:));
        bc_f = squeeze(all_bc(:,f,:));
        par_f = squeeze(all_par(:,f,:));
        ec_f = squeeze(all_ec(:,f,:));
        clust_f = squeeze(all_clust(:,f,:));
        cc_reg_f = squeeze(cc_reg(:,f,:));
       % le_f = squeeze(all_le(:,f,:));
        
        % Loop over each permutation
        for i_p = 1:n_perm
            c_c_f_p = squeeze(c_c_f(:,i_p));
            ns_f_p = squeeze(ns_f(:,i_p));
            bc_f_p = squeeze(bc_f(:,i_p));
            par_f_p = squeeze(par_f(:,i_p));
            ec_f_p = squeeze(ec_f(:,i_p));
            clust_f_p = squeeze(clust_f(:,i_p));
            cc_reg_f_p = squeeze(cc_reg_f(:,i_p));
           % le_f_p = squeeze(le_f(:,i_p));

            %% Find most synchronizing node

            % Get the identity of the most synchronizing node (the one we would
            % tell the surgeons to resect)
            [~,ch_most_sync] = min(c_c_f_p);
            most_sync(f,i_p) = ch_most_sync;
            

            %% Do Spearman rank for nodal measures
            [rho_cc(f,i_p),~] = doStats(c_c,c_c_f_p);
            [rho_ns(f,i_p),~] = doStats(ns,ns_f_p);
            [rho_bc(f,i_p),~] = doStats(bc,bc_f_p);
            [rho_par(f,i_p),~] = doStats(par,par_f_p);
            [rho_ec(f,i_p),~] = doStats(ec,ec_f_p);
            [rho_clust(f,i_p),~] = doStats(clust,clust_f_p);
            [rho_cc_reg(f,i_p),~] = doStats(cc_regional,cc_reg_f_p);
           % [rho_le(f,i_p),~] = doStats(le,le_f_p);
            
            %% Get rho just for electrodes in the resection zone
            if isempty(pt(whichPt).resec) == 1
                rho_cc_resec(f,i_p) = nan;
                rho_bc_resec(f,i_p) = nan;
                rho_ns_resec(f,i_p) = nan;
                rho_par_resec(f,i_p) = nan;
                rho_ec_resec(f,i_p) = nan;
                rho_clust_resec(f,i_p) = nan;
                %rho_le_resec(f,i_p) = nan;
            else
                resec = pt(whichPt).resec.nums;
                [rho_ns_resec(f,i_p),~] = doStats(ns(resec),ns_f_p(resec));
                [rho_bc_resec(f,i_p),~] = doStats(bc(resec),bc_f_p(resec));
                [rho_cc_resec(f,i_p),~] = doStats(c_c(resec),c_c_f_p(resec));
                [rho_par_resec(f,i_p),~] = doStats(par(resec),par_f_p(resec));
                [rho_ec_resec(f,i_p),~] = doStats(ec(resec),ec_f_p(resec));
                [rho_clust_resec(f,i_p),~] = doStats(clust(resec),clust_f_p(resec));
               % [rho_le_resec(f,i_p),~] = doStats(le(resec),le_f_p(resec));
            end
                

        end

    end

    %% Agreement measures, which I will include in both analyses
    % Global measures of agreement: relative difference
    rel_sync = (all_sync-sync)/sync;
    rel_eff = (all_eff-eff)/eff;
    rel_trans = (all_trans-trans)/trans;
    rel_sync_norm = (all_sync_norm-sync_norm)/sync_norm;
    rel_eff_norm = (all_eff_norm-eff_norm)/eff_norm;
    rel_trans_norm = (all_trans_norm-trans_norm)/trans_norm;
    

    % Nodal measures: "average" the SRCs. To average the SRCs, I
    % apply a Fisher's transformation, average the z values, and then back
    % transform to get an average rho.
    rho_mean_cc = average_rho(rho_cc,2);
    rho_mean_ns = average_rho(rho_ns,2);
    rho_mean_bc = average_rho(rho_bc,2);
    rho_mean_par = average_rho(rho_par,2);
    rho_mean_ec = average_rho(rho_ec,2);
    rho_mean_clust = average_rho(rho_clust,2);
    rho_mean_cc_reg = average_rho(rho_cc_reg,2);
   % rho_mean_le = average_rho(rho_le,2);
     
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
        par_rel_std = rel_std_nodal(all_par,par);
        ec_rel_std = rel_std_nodal(all_ec,ec);
        clust_rel_std = rel_std_nodal(all_clust,clust);
        cc_reg_rel_std = rel_std_nodal(cc_reg,cc_regional);
       % le_rel_std = rel_std_nodal(all_le,le);
        
        %% Global measures
        % Global measures: we will do relative std (std across
        % permutations divided by std across patients). And so we will
        % record the standard deviation for now.
        
       
        %% Fill up structures
        stats(whichPt).name = name;
 
        % control centrality
        stats(whichPt).(contig_text).(sec_text).cc.true = c_c;
        stats(whichPt).(contig_text).(sec_text).cc.rel_std = cc_rel_std;
        stats(whichPt).(contig_text).(sec_text).cc.rho_mean = rho_mean_cc;
        
        % regional control centrality
        stats(whichPt).(contig_text).(sec_text).cc_reg.true = cc_regional;
        stats(whichPt).(contig_text).(sec_text).cc_reg.rel_std = cc_reg_rel_std;
        stats(whichPt).(contig_text).(sec_text).cc_reg.rho_mean = rho_mean_cc_reg;
        
        % Most synchronizing electrode
        stats(whichPt).(contig_text).(sec_text).min_cc.true = min_cc_true;
        
        % Electrodes in most synchronizing contiguous region
        stats(whichPt).(contig_text).(sec_text).regional_cc.true = elecs_regional_min;
        
        % Electrodes in various percentiles of most synchronizing electrode
        % and region across resampling
        for i = [70 80 90 95]
            elecs_single = get_perc_elecs(most_sync(4,:),i);
            single_text = sprintf('single_%d',i);
            reg_text = sprintf('regional_%d',i);
            if isempty(elecs_min) == 0
                elecs_regional = get_perc_elecs(elecs_min,i);
            else
                elecs_regional = nan;
            end
            stats(whichPt).(contig_text).(sec_text).min_cc_elecs.(single_text) = elecs_single;
            stats(whichPt).(contig_text).(sec_text).min_cc_elecs.(reg_text) = elecs_regional;
        end

        % node strength
        stats(whichPt).(contig_text).(sec_text).ns.rel_std = ns_rel_std;
        stats(whichPt).(contig_text).(sec_text).ns.rho_mean = rho_mean_ns;

        % betweenness centrality
        stats(whichPt).(contig_text).(sec_text).bc.rel_std = bc_rel_std;
        stats(whichPt).(contig_text).(sec_text).bc.rho_mean = rho_mean_bc;
        
        % Participation coeff
        stats(whichPt).(contig_text).(sec_text).par.rel_std = par_rel_std;
        stats(whichPt).(contig_text).(sec_text).par.rho_mean = rho_mean_par;
        
        % Eigenvector centrality
        stats(whichPt).(contig_text).(sec_text).ec.rel_std = ec_rel_std;
        stats(whichPt).(contig_text).(sec_text).ec.rho_mean = rho_mean_ec;
        
        % Local efficiency
       % stats(whichPt).(contig_text).(sec_text).le.rel_std = le_rel_std;
       % stats(whichPt).(contig_text).(sec_text).le.rho_mean = rho_mean_le;
        
        % clustering coefficient
        stats(whichPt).(contig_text).(sec_text).clust.rel_std = clust_rel_std;
        stats(whichPt).(contig_text).(sec_text).clust.rho_mean = rho_mean_clust;

        % synchronizability
        stats(whichPt).(contig_text).(sec_text).sync.std = std(all_sync,0,2);
        stats(whichPt).(contig_text).(sec_text).sync.true = sync;
        stats(whichPt).(contig_text).(sec_text).sync.rel_diff = rel_sync;
        stats(whichPt).(contig_text).(sec_text).sync.rel_diff_norm = rel_sync_norm;
        
        % transitivity
        stats(whichPt).(contig_text).(sec_text).trans.std = std(all_trans,0,2);
        stats(whichPt).(contig_text).(sec_text).trans.true = trans;
        stats(whichPt).(contig_text).(sec_text).trans.rel_diff = rel_trans;
        stats(whichPt).(contig_text).(sec_text).trans.rel_diff_norm = rel_trans_norm;

        % efficiency
        stats(whichPt).(contig_text).(sec_text).eff.std = std(all_eff,0,2);
        stats(whichPt).(contig_text).(sec_text).eff.true = eff;
        stats(whichPt).(contig_text).(sec_text).eff.rel_diff = rel_eff;
        stats(whichPt).(contig_text).(sec_text).eff.rel_diff_norm = rel_eff_norm;

        if doSave == 1
            save([resultsFolder,'basic_metrics/stats.mat'],'stats');
        end
    elseif do_soz_analysis == 1
        %% Do the analysis of dependence of agreement on distance from important things
       
        % Nodal
        soz(whichPt).(contig_text).(sec_text).rho_cc = rho_cc;
        soz(whichPt).(contig_text).(sec_text).rho_bc = rho_bc;
        soz(whichPt).(contig_text).(sec_text).rho_ns = rho_ns;
        soz(whichPt).(contig_text).(sec_text).rho_ec = rho_ec;
        soz(whichPt).(contig_text).(sec_text).rho_par = rho_par;
        soz(whichPt).(contig_text).(sec_text).rho_clust = rho_clust;
        soz(whichPt).(contig_text).(sec_text).rho_cc_reg = rho_cc_reg;
       % soz(whichPt).(contig_text).(sec_text).rho_le = rho_le;
        
        % Global
        soz(whichPt).(contig_text).(sec_text).sync = rel_sync;
        soz(whichPt).(contig_text).(sec_text).eff = rel_eff;
        soz(whichPt).(contig_text).(sec_text).trans = rel_trans;
        soz(whichPt).(contig_text).(sec_text).dist_soz = dist_soz;
        soz(whichPt).(contig_text).(sec_text).dist_resec = dist_resec;
        soz(whichPt).(contig_text).(sec_text).overlap_soz = overlap_soz;
        soz(whichPt).(contig_text).(sec_text).overlap_resec = overlap_resec;
        soz(whichPt).(contig_text).(sec_text).par_removed = avg_par_removed;
        soz(whichPt).(contig_text).(sec_text).bc_removed = avg_bc_removed;
        
        soz(whichPt).(contig_text).(sec_text).rho_cc_resec = rho_cc_resec;
        soz(whichPt).(contig_text).(sec_text).rho_bc_resec = rho_bc_resec;
        soz(whichPt).(contig_text).(sec_text).rho_ns_resec = rho_ns_resec;
        soz(whichPt).(contig_text).(sec_text).rho_par_resec = rho_par_resec;
        soz(whichPt).(contig_text).(sec_text).rho_ec_resec = rho_ec_resec;
        soz(whichPt).(contig_text).(sec_text).rho_clust_resec = rho_clust_resec;
       % soz(whichPt).(contig_text).(sec_text).rho_le_resec = rho_le_resec;
        if doSave == 1
            save([resultsFolder,'basic_metrics/soz.mat'],'soz');
        end
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
        
    % Std across electrodes of the true metric
    std_across_ch = nanstd(true_metric,0,1);

    % Relative std
    rel_std = std_across_perm./std_across_ch;

end