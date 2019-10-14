function network_stats(whichPts,do_soz_analysis,which_sz)

%{
This function takes an adjacency matrix A, calculates global and nodal
network measures, and then randomly resamples the network, retaining
specified fractions of electrodes, and then gets statistics on how the
network measures change.


for whichPts, should put in 1:33 to do the whole run
%}

total_time = 0;
rng('default')


%% Parameters

% Which sz: if 1, the first the patient has, if 2, the second

% Save the output? (Should be yes)
doSave = 1;

% add to existing stats array? Or re-write with new stats array? Usually
% yes
merge = 1;

% do_soz_analysis: 1 if doing SOZ analysis, 0 if doing main analysis, 2 if
% doing new soz overlap analysis
% e_f: What fraction of nodes to retain
% contigs: 1 means contiguous set of electrodes, 0 means random electrodes
if do_soz_analysis == 1
    e_f = [0.2 0.4 0.6 0.8 1]; % JUST CHANGED THIS primary way, systematically remove electrode and its N
    %nearest neighbors, amounting to 20% of total number of electrodes
   % e_f = 999; % other way, we are just removing single channel
    contigs = 1;
elseif do_soz_analysis == 0
    e_f = [0.2 0.4 0.6 0.8 1];
    contigs = [0 1];
elseif do_soz_analysis == 2 % compare targeting soz to sparing soz
    contigs = [2 3];
    e_f = nan;
elseif do_soz_analysis == 3 % compare targeting soz to random
    contigs = [4 3];
    e_f = nan;
end
n_f = length(e_f);



%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);

if which_sz == 1
    extra = '';
else
    extra = '2';
end


%% Load the output structure to add more info to it
if merge == 1
    if do_soz_analysis == 1
        if exist([resultsFolder,'basic_metrics/soz',extra,'.mat'],'file') ~= 0
            load([resultsFolder,'basic_metrics/soz',extra,'.mat']);
        else
            soz = struct;
        end
    elseif do_soz_analysis == 0
        if exist([resultsFolder,'basic_metrics/stats',extra,'.mat'],'file') ~= 0
            load([resultsFolder,'basic_metrics/stats',extra,'.mat']);
        else
            stats = struct;
        end
    elseif do_soz_analysis == 2
        if exist([resultsFolder,'basic_metrics/soz_overlap',extra,'.mat'],'file') ~= 0
            fprintf('Found existing file, loading...\n');
            load([resultsFolder,'basic_metrics/soz_overlap',extra,'.mat']);
        else
            soz_overlap = struct;
        end
    elseif do_soz_analysis == 3
        if exist([resultsFolder,'basic_metrics/soz_overlap_random',extra,'.mat'],'file') ~= 0
            fprintf('Found existing file, loading...\n');
            load([resultsFolder,'basic_metrics/soz_overlap_random',extra,'.mat']);
        else
            soz_overlap = struct;
        end
    end
else
    if do_soz_analysis == 1
        soz = struct;
    elseif do_soz_analysis == 0
        stats = struct;
    elseif do_soz_analysis == 2 || do_soz_analysis == 3
        soz_overlap = struct;
    end
end

%% which frequencies to do
freq_cell = {'high_gamma','beta'};

%% Loop through patients, times, frequencies, and whether contig or random electrodes
for ff = 1:length(freq_cell)
    freq = freq_cell{ff};
    fprintf('Doing %s\n',freq);
for which_sec = [0 -10 -5 5 10] % 0 means EEC, -5 is 5 seconds before
    fprintf('Doing %d second\n',which_sec);

    
% Loop through patients
for whichPt = whichPts

for contig = contigs % random or contiguous electrodes
    fprintf('Doing contig %d\n',contig);
    tic
    
    % Skip if all electrode locations are -1 (means we don't have electrode
    % locations)
    if unique(pt(whichPt).new_elecs.locs) == -1
        continue
    end
    
    if contig == 1
        % Here, not taking random samples, but rather systematically going
        % through each electrode and its N nearest neighbors
        n_perm = length(pt(whichPt).new_elecs.electrodes);
    elseif contig == 0 || contig == 2
        % Take 1000 random permutations
        n_perm = 1e3;
    elseif contig == 3
        n_perm = 1;
    end

    % Make result folder
    name = pt(whichPt).name;
    
    fprintf('Doing %s\n',name);
    
    outFolder = [resultsFolder,'basic_metrics/',name,'/'];
    %{
    if exist(outFolder,'dir') == 0
        mkdir(outFolder);
    end
    %}

    if contig == 1
        contig_text = 'contiguous';
    elseif contig == 0
        contig_text = 'random';
    elseif contig == 2
        contig_text = 'not_soz';
    elseif contig == 3
        contig_text = 'soz';
    elseif contig == 4
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
        elseif do_soz_analysis == 2 || do_soz_analysis == 3
            stats = soz_overlap;
        end
        if length(stats) >= whichPt
            if isfield(stats(whichPt),(freq)) == 1
                if isfield(stats(whichPt).(freq),(contig_text)) == 1
                    if isfield(stats(whichPt).(freq).(contig_text),(sec_text)) == 1

                        fprintf('Did %s, skipping\n',name);
                        continue
                    end
                end 
            end
        end
    end

    %% Get adjacency matrix
    [adj,~,sz_num] = reconcileAdj(pt,whichPt,which_sz);
    if isempty(adj) == 1
        fprintf('Cannot do %s\n\n',name);
        continue;
    end

    %% Get appropriate frequency band
    if strcmp(freq,'high_gamma') == 1
        A_all = adj(4).data;
        if contains(adj(4).name,'highgamma') == 0
            error('This isn''t gamma!'\n');
        end
    elseif strcmp(freq,'beta') == 1
        A_all = adj(2).data;
        if contains(adj(2).name,'beta') == 0
            error('This isn''t beta!\n');
        end
    end

    %% Get which time
    % Start with the middle and add which second
    if ceil(size(A_all,1)/2)+which_sec <= 0, continue; end
    if ceil(size(A_all,1)/2)+which_sec > size(A_all,1), continue; end
    A = squeeze(A_all(ceil(size(A_all,1)/2)+which_sec,:,:));
    if sum(sum(isnan(A))) == sum(sum(ones(size(A))))
        continue
    end

    %% Get true metrics
    fprintf('Getting true metrics\n');
    
    %% Control centrality
    c_c = control_centrality(A);
    
    cc_fake = nan(100,1);
    for i = 1:100
        cc_fake(i) = mean(control_centrality(generate_fake_graph(A)));
    end
    m_cc_fake = mean(control_centrality(generate_fake_graph(A)));
    cc_norm = c_c./m_cc_fake;
    
    % Get identity of node with lowest control centrality
    [~,min_cc_true] = min(c_c);
    
    % Get locations
    locs = pt(whichPt).new_elecs.locs; % all electrode locations
    
    %% Second order control centrality
    %{
    socc = second_order_cc(A);
    
    figure
    subplot(2,1,1)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
    hold on
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,c_c,'filled');
    set(gca,'clim',prctile(c_c,[10 90]));
    
    subplot(2,1,2)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
    hold on
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,1./socc,'filled');
    set(gca,'clim',prctile(1./socc,[10 90]));
    %}
    
    %% Regional control centrality
    % the change in synchronizability when you remove an electrode and its
    % N nearest neighbors, where N = the number of resected electrodes
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
    
    % Control for graph size by dividing by an average fake graph of the
    % same size, obtained by permuting entries of the upper triangular of
    % the true A, and reflecting over to get an undirected graph.
    sync_fake = nan(100,1);
    for i = 1:100
        sync_fake(i) = synchronizability(generate_fake_graph(A));
    end
    m_sync_fake = mean(sync_fake);
    sync_norm = sync/m_sync_fake;

    %% Get true betweenness centrality
    bc = betweenness_centrality(A,1);
    bc_fake = nan(100,1);
    for i = 1:100
        bc_fake(i) = mean(betweenness_centrality(generate_fake_graph(A),1));
    end
    m_bc_fake = mean(bc_fake);
    bc_norm = bc./m_bc_fake;
    
    % Get identity of node with highest betweenness centrality
    [~,max_bc_true] = max(bc);
    
    %% Get true eigenvector centrality
    ec = eigenvector_centrality_und(A);
    ec_fake = nan(100,1);
    for i = 1:100
        ec_fake(i) = mean(eigenvector_centrality_und(generate_fake_graph(A)));
    end
    m_ec_fake = mean(ec_fake);
    ec_norm = ec./m_ec_fake;
    
    % Get identity of node with highest eigenvector centrality
    [~,max_ec_true] = max(ec);
    
    %% Get true clustering coefficient
    clust = clustering_coef_wu(A);
    clust_fake = nan(100,1);
    for i = 1:100
        clust_fake(i) = mean(clustering_coef_wu(generate_fake_graph(A)));
    end
    m_clust_fake = mean(clust_fake);
    clust_norm = clust./m_clust_fake;
    
    % Get identity of node with highest clustering coefficient
    [~,max_clust_true] = max(clust);
    
    %% Get local efficiency
   % le = efficiency_wei(A,1);

    %% Get true node strength
    ns = node_strength(A);
    ns_fake = nan(100,1);
    for i = 1:100
        ns_fake(i) = mean(node_strength(generate_fake_graph(A)));
    end
    m_ns_fake = mean(ns_fake);
    ns_norm = ns./m_ns_fake;
    
    % Get identity of node with highest node strength
    [~,max_ns_true] = max(ns);
    
    %% Get true participation coefficient
    [Ci,~] = modularity_und(A);
    par = participation_coef(A,Ci);

    %% Get true global efficiency and efficiency for fake network
    eff = efficiency_wei(A, 0);
    eff_fake = nan(100,1);
    for i = 1:100
        eff_fake(i) = efficiency_wei(generate_fake_graph(A),0);
    end
    m_eff_fake = mean(eff_fake);
    eff_norm = eff/m_eff_fake;
    
    %% Get true transitivity and normalized transitivity
    trans = transitivity_wu(A);
    trans_fake = nan(100,1);
    for i = 1:100
        trans_fake(i) = transitivity_wu(generate_fake_graph(A));
    end
    m_trans_fake = mean(trans_fake);
    trans_norm = trans/m_trans_fake;

    
    fprintf('Got true metrics, now resampling network...\n');
    
    %% Resample network and get new metrics
    % all_c_c is nch x n_f x n_perm size matrix
    [all_c_c,all_ns,all_bc,all_sync,all_eff,overlap_soz,dist_soz,...
        overlap_resec,dist_resec,elecs_min,...
        all_par,all_trans,avg_par_removed,avg_bc_removed,...
        all_sync_norm,all_eff_norm,all_trans_norm,all_ec,...
        all_clust,all_le,cc_reg,dist_nearest_resec,sz_soz_dist,...
        all_cc_norm,all_ns_norm,all_bc_norm,all_ec_norm,all_clust_norm] = ...
        resampleNetwork(A,n_perm,e_f,contig,pt,whichPt,adj,sz_num);

    %% Initialize arrays to compare old to new metrics

    % Control centrality stuff
    rho_cc = zeros(n_f,n_perm);
    most_sync = zeros(n_f,n_perm);

    % Betweenness centrality stuff
    rho_bc = zeros(n_f,n_perm);
    most_bc = zeros(n_f,n_perm);

    % Node strength centrality stuff
    rho_ns = zeros(n_f,n_perm);
    most_ns = zeros(n_f,n_perm);
    
    % Participation coefficient stuff
    rho_par = zeros(n_f,n_perm);
    
    % Eigenvalue centality stuff
    rho_ec = zeros(n_f,n_perm);
    most_ec = zeros(n_f,n_perm);
    
    % Local efficiency stuff
    rho_le = zeros(n_f,n_perm);
    
    % Clustering coefficient stuff
    rho_clust = zeros(n_f,n_perm);
    most_clust = zeros(n_f,n_perm);
    
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
    
    % Changing identity of the hub electrode
    same_most_ns = zeros(n_f,n_perm);
    same_most_clust = zeros(n_f,n_perm);
    same_most_ec = zeros(n_f,n_perm);
    same_most_bc = zeros(n_f,n_perm);
    same_most_sync = zeros(n_f,n_perm);
    
    %% Loop over each fraction and get various stats
    for f = 1:n_f

        % The metrics for that number of removed electrodes
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
            
            % The metrics for that specific permutation
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
            
            [~,most_ns(f,i_p)] = max(ns_f_p);
            [~,most_bc(f,i_p)] = max(bc_f_p);
            [~,most_ec(f,i_p)] = max(ec_f_p);
            [~,most_clust(f,i_p)] = max(clust_f_p);
            
            % Decide if this hub is the same as in the original network
            if most_ns(f,i_p) == max_ns_true
                same_most_ns(f,i_p) = 1;
            end
            
            if most_bc(f,i_p) == max_bc_true
                same_most_bc(f,i_p) = 1;
            end
            
            if most_ec(f,i_p) == max_ec_true
                same_most_ec(f,i_p) = 1;
            end
                
            if most_clust(f,i_p) == max_clust_true
                same_most_clust(f,i_p) = 1;
            end
            
            if most_sync(f,i_p) == min_cc_true
                same_most_sync(f,i_p) = 1;
            end
            

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
        
        %% Get alternate norm
        other_cc_norm = all_c_c./squeeze(nanmean(nanmean(all_c_c,3),1));
        other_bc_norm = all_bc./squeeze(nanmean(nanmean(all_bc,3),1));
        other_ns_norm = all_ns./squeeze(nanmean(nanmean(all_ns,3),1));
        other_ec_norm = all_ec./squeeze(nanmean(nanmean(all_ec,3),1));
        other_clust_norm = all_clust./squeeze(nanmean(nanmean(all_clust,3),1));
       
        %% Nodal measure of variability
        % Nodal measure: relative std (std across permuations divided by
        % std across electrodes)
        %{
        cc_rel_std = rel_std_nodal(all_c_c,c_c);
        ns_rel_std = rel_std_nodal(all_ns,ns);
        bc_rel_std = rel_std_nodal(all_bc,bc);
        par_rel_std = rel_std_nodal(all_par,par);
        ec_rel_std = rel_std_nodal(all_ec,ec);
        clust_rel_std = rel_std_nodal(all_clust,clust);
        cc_reg_rel_std = rel_std_nodal(cc_reg,cc_regional);
       % le_rel_std = rel_std_nodal(all_le,le);
        %}
        
       
        
        % Reliability
        cc_rel = reliability_nodal(all_c_c,c_c);
        ns_rel = reliability_nodal(all_ns,ns);
        bc_rel = reliability_nodal(all_bc,bc);
        par_rel = reliability_nodal(all_par,par);
        ec_rel = reliability_nodal(all_ec,ec);
        clust_rel = reliability_nodal(all_clust,clust);
        cc_reg_rel = reliability_nodal(cc_reg,cc_regional);
        
        cc_rel_norm = reliability_nodal(all_cc_norm,cc_norm);
        ns_rel_norm = reliability_nodal(all_ns_norm,ns_norm);
        bc_rel_norm = reliability_nodal(all_bc_norm,bc_norm);
        ec_rel_norm = reliability_nodal(all_ec_norm,ec_norm);
        clust_rel_norm = reliability_nodal(all_clust_norm,clust_norm);
        
        other_cc_rel_norm = reliability_nodal(other_cc_norm,cc_norm);
        other_ns_rel_norm = reliability_nodal(other_ns_norm,ns_norm);
        other_bc_rel_norm = reliability_nodal(other_bc_norm,bc_norm);
        other_ec_rel_norm = reliability_nodal(other_ec_norm,ec_norm);
        other_clust_rel_norm = reliability_nodal(other_clust_norm,clust_norm);
        
        
        % Alternate reliability
        cc_rel_alt = alt_rel_nodal(all_c_c);
        ns_rel_alt = alt_rel_nodal(all_ns);
        bc_rel_alt = alt_rel_nodal(all_bc);
        par_rel_alt = alt_rel_nodal(all_par);
        ec_rel_alt = alt_rel_nodal(all_ec);
        clust_rel_alt = alt_rel_nodal(all_clust);
        cc_reg_rel_alt = alt_rel_nodal(cc_reg);
        
        %% Global measures
        % Global measures: we will do global reliability so for now store
        % the std across permutations
        
       
        %% Fill up structures
        stats(whichPt).name = name;
         
        % control centrality
        stats(whichPt).(freq).(contig_text).(sec_text).cc.true = c_c;
        stats(whichPt).(freq).(contig_text).(sec_text).cc.rel = cc_rel;
        stats(whichPt).(freq).(contig_text).(sec_text).cc.rel_alt = cc_rel_alt;
        stats(whichPt).(freq).(contig_text).(sec_text).cc.rel_norm = cc_rel_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).cc.other_rel_norm = other_cc_rel_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).cc.rho_mean = rho_mean_cc;
        stats(whichPt).(freq).(contig_text).(sec_text).cc.same_hub = same_most_sync;
        
        % regional control centrality
        stats(whichPt).(freq).(contig_text).(sec_text).cc_reg.true = cc_regional;
        stats(whichPt).(freq).(contig_text).(sec_text).cc_reg.rel = cc_reg_rel;
        stats(whichPt).(freq).(contig_text).(sec_text).cc_reg.rel_alt = cc_reg_rel_alt;
        stats(whichPt).(freq).(contig_text).(sec_text).cc_reg.rho_mean = rho_mean_cc_reg;
        
        % Most synchronizing electrode
        stats(whichPt).(freq).(contig_text).(sec_text).min_cc.true = min_cc_true;
        
        % Electrodes in most synchronizing contiguous region
        stats(whichPt).(freq).(contig_text).(sec_text).regional_cc.true = elecs_regional_min;
        
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
            stats(whichPt).(freq).(contig_text).(sec_text).min_cc_elecs.(single_text) = elecs_single;
            stats(whichPt).(freq).(contig_text).(sec_text).min_cc_elecs.(reg_text) = elecs_regional;
            
            stats(whichPt).(freq).(contig_text).(sec_text).cc.(single_text) = get_perc_elecs(most_sync(4,:),i);
            stats(whichPt).(freq).(contig_text).(sec_text).ns.(single_text) = get_perc_elecs(most_ns(4,:),i);
            stats(whichPt).(freq).(contig_text).(sec_text).bc.(single_text) = get_perc_elecs(most_bc(4,:),i);
            stats(whichPt).(freq).(contig_text).(sec_text).ec.(single_text) = get_perc_elecs(most_ec(4,:),i);
            stats(whichPt).(freq).(contig_text).(sec_text).clust.(single_text) = get_perc_elecs(most_clust(4,:),i);
        end
        
        [~,stats(whichPt).(freq).(contig_text).(sec_text).ns.true] = max(ns);
        [~,stats(whichPt).(freq).(contig_text).(sec_text).bc.true] = max(bc);
        [~,stats(whichPt).(freq).(contig_text).(sec_text).ec.true] = max(ec);
        [~,stats(whichPt).(freq).(contig_text).(sec_text).clust.true] = max(clust);


        % node strength
        
        stats(whichPt).(freq).(contig_text).(sec_text).ns.rel = ns_rel;
        stats(whichPt).(freq).(contig_text).(sec_text).ns.rel_alt = ns_rel_alt;
        stats(whichPt).(freq).(contig_text).(sec_text).ns.rel_norm = ns_rel_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).ns.other_rel_norm = other_ns_rel_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).ns.rho_mean = rho_mean_ns;
        stats(whichPt).(freq).(contig_text).(sec_text).ns.same_hub = same_most_ns;

        % betweenness centrality
        stats(whichPt).(freq).(contig_text).(sec_text).bc.rel = bc_rel;
        stats(whichPt).(freq).(contig_text).(sec_text).bc.rel_alt = bc_rel_alt;
        stats(whichPt).(freq).(contig_text).(sec_text).bc.rel_norm = bc_rel_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).bc.other_rel_norm = other_bc_rel_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).bc.rho_mean = rho_mean_bc;
        stats(whichPt).(freq).(contig_text).(sec_text).bc.same_hub = same_most_bc;
        
        % Participation coeff
        stats(whichPt).(freq).(contig_text).(sec_text).par.rel = par_rel;
        stats(whichPt).(freq).(contig_text).(sec_text).par.rel_alt = par_rel_alt;
        stats(whichPt).(freq).(contig_text).(sec_text).par.rho_mean = rho_mean_par;
        
        % Eigenvector centrality
        stats(whichPt).(freq).(contig_text).(sec_text).ec.rel = ec_rel;
        stats(whichPt).(freq).(contig_text).(sec_text).ec.rel_alt = ec_rel_alt;
        stats(whichPt).(freq).(contig_text).(sec_text).ec.rel_norm = ec_rel_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).ec.other_rel_norm = other_ec_rel_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).ec.rho_mean = rho_mean_ec;
        stats(whichPt).(freq).(contig_text).(sec_text).ec.same_hub = same_most_ec;
        
        % Local efficiency
       % stats(whichPt).(contig_text).(sec_text).le.rel_std = le_rel_std;
       % stats(whichPt).(contig_text).(sec_text).le.rho_mean = rho_mean_le;
        
        % clustering coefficient
        stats(whichPt).(freq).(contig_text).(sec_text).clust.rel = clust_rel;
        stats(whichPt).(freq).(contig_text).(sec_text).clust.rel_alt = clust_rel_alt;
        stats(whichPt).(freq).(contig_text).(sec_text).clust.rel_norm = clust_rel_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).clust.other_rel_norm = other_clust_rel_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).clust.rho_mean = rho_mean_clust;
        stats(whichPt).(freq).(contig_text).(sec_text).clust.same_hub = same_most_clust;

        % synchronizability
        stats(whichPt).(freq).(contig_text).(sec_text).sync.std = std(all_sync,0,2);
        stats(whichPt).(freq).(contig_text).(sec_text).sync.std_norm = std(all_sync_norm,0,2);
        stats(whichPt).(freq).(contig_text).(sec_text).sync.true = sync;
        stats(whichPt).(freq).(contig_text).(sec_text).sync.true_norm = sync_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).sync.rel_diff = rel_sync;
        stats(whichPt).(freq).(contig_text).(sec_text).sync.rel_diff_norm = rel_sync_norm;
        
        % transitivity
        stats(whichPt).(freq).(contig_text).(sec_text).trans.std = std(all_trans,0,2);
        stats(whichPt).(freq).(contig_text).(sec_text).trans.std_norm = std(all_trans_norm,0,2);
        stats(whichPt).(freq).(contig_text).(sec_text).trans.true = trans;
        stats(whichPt).(freq).(contig_text).(sec_text).trans.true_norm = trans_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).trans.rel_diff = rel_trans;
        stats(whichPt).(freq).(contig_text).(sec_text).trans.rel_diff_norm = rel_trans_norm;

        % efficiency
        stats(whichPt).(freq).(contig_text).(sec_text).eff.std = std(all_eff,0,2);
        stats(whichPt).(freq).(contig_text).(sec_text).eff.std_norm = std(all_eff_norm,0,2);
        stats(whichPt).(freq).(contig_text).(sec_text).eff.true = eff;
        stats(whichPt).(freq).(contig_text).(sec_text).eff.true_norm = eff_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).eff.rel_diff = rel_eff;
        stats(whichPt).(freq).(contig_text).(sec_text).eff.rel_diff_norm = rel_eff_norm;
        
        
        % Save all global measures
        stats(whichPt).(freq).(contig_text).(sec_text).sync.all = all_sync;
        stats(whichPt).(freq).(contig_text).(sec_text).eff.all = all_eff;
        stats(whichPt).(freq).(contig_text).(sec_text).trans.all = all_trans;
        
        stats(whichPt).(freq).(contig_text).(sec_text).sync.all_norm = all_sync_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).eff.all_norm = all_eff_norm;
        stats(whichPt).(freq).(contig_text).(sec_text).trans.all_norm = all_trans_norm;

        if doSave == 1
            save([resultsFolder,'basic_metrics/stats',extra,'.mat'],'stats');
        end
    elseif do_soz_analysis == 1 || do_soz_analysis == 2 || do_soz_analysis == 3
        %% Do the analysis of dependence of agreement on distance from important things
        if do_soz_analysis == 2 || do_soz_analysis == 3
            soz = soz_overlap;
        end
        
        % Change hub
        soz(whichPt).(freq).(contig_text).(sec_text).same_hub.ec = same_most_ec;
        soz(whichPt).(freq).(contig_text).(sec_text).same_hub.bc = same_most_bc;
        soz(whichPt).(freq).(contig_text).(sec_text).same_hub.ns = same_most_ns;
        soz(whichPt).(freq).(contig_text).(sec_text).same_hub.clust = same_most_clust;
        soz(whichPt).(freq).(contig_text).(sec_text).same_hub.cc = same_most_sync;
        
        % Nodal
        soz(whichPt).(freq).(contig_text).(sec_text).rho_cc = rho_cc;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_bc = rho_bc;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_ns = rho_ns;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_ec = rho_ec;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_par = rho_par;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_clust = rho_clust;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_cc_reg = rho_cc_reg;
       % soz(whichPt).(contig_text).(sec_text).rho_le = rho_le;
        
        % Global
        soz(whichPt).(freq).(contig_text).(sec_text).sync = rel_sync;
        soz(whichPt).(freq).(contig_text).(sec_text).eff = rel_eff;
        soz(whichPt).(freq).(contig_text).(sec_text).trans = rel_trans;
        
        % Distance from/overlap with SOZ/resection zone
        soz(whichPt).(freq).(contig_text).(sec_text).sz_soz_dist = sz_soz_dist;
        soz(whichPt).(freq).(contig_text).(sec_text).dist_soz = dist_soz;
        soz(whichPt).(freq).(contig_text).(sec_text).dist_resec = dist_resec;
        soz(whichPt).(freq).(contig_text).(sec_text).overlap_soz = overlap_soz;
        soz(whichPt).(freq).(contig_text).(sec_text).overlap_resec = overlap_resec;
        soz(whichPt).(freq).(contig_text).(sec_text).par_removed = avg_par_removed;
        soz(whichPt).(freq).(contig_text).(sec_text).bc_removed = avg_bc_removed;
        soz(whichPt).(freq).(contig_text).(sec_text).dist_nearest_resec = dist_nearest_resec;
        
        % Agreement of nodal measures for electrodes in resection zone
        soz(whichPt).(freq).(contig_text).(sec_text).rho_cc_resec = rho_cc_resec;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_bc_resec = rho_bc_resec;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_ns_resec = rho_ns_resec;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_par_resec = rho_par_resec;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_ec_resec = rho_ec_resec;
        soz(whichPt).(freq).(contig_text).(sec_text).rho_clust_resec = rho_clust_resec;
       % soz(whichPt).(contig_text).(sec_text).rho_le_resec = rho_le_resec;
        if doSave == 1 
            if do_soz_analysis == 1
                save([resultsFolder,'basic_metrics/soz',extra,'.mat'],'soz');
            elseif do_soz_analysis == 2
                soz_overlap = soz;
                save([resultsFolder,'basic_metrics/soz_overlap',extra,'.mat'],'soz_overlap');
            elseif do_soz_analysis == 3
                soz_overlap = soz;
                save([resultsFolder,'basic_metrics/soz_overlap_random',extra,'.mat'],'soz_overlap');
            end
        end
    end

    %% Plots
    % Do patient-specific plots (usually no)
    doPlots = 0;
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
    
   t = toc;
   fprintf('Elapsed time for %s was %1.2f minutes.\n\n',name,t/60); 
   total_time =  total_time + t;
    
end

end
end
end

fprintf('The total time was %1.2f hours.\n\n',total_time/3600);

end
