function [all_c_c,all_ns,all_bc,all_sync,all_eff,overlap_soz,dist_soz,...
    overlap_resec,dist_resec,elecs_min,all_par,all_trans,...
    avg_par_removed,avg_bc_removed,all_sync_norm,all_eff_norm,all_trans_norm,...
    all_ec,all_clust,all_le,cc_reg,dist_nearest_resec,sz_soz_dist,...
    all_cc_norm,all_ns_norm,all_bc_norm,all_ec_norm,all_clust_norm] = ...
    resampleNetwork(A,n_perm,e_f,contig,pt,whichPt,adj,sz_num)

%{
This function resamples the network by removing a fraction of nodes and then
recalculates various network metrics    
%}
    
    
%% Get basic numbers

doPlot = 0; % Plot overlap of electrodes?

% Get electrode locs
locs = pt(whichPt).new_elecs.locs;

% Get soz electrodes
soz = pt(whichPt).soz.nums;
sz_soz = pt(whichPt).sz(sz_num).nums;

% Get resected electrodes
if isempty(pt(whichPt).resec) == 0
    resec = pt(whichPt).resec.nums;
else
    resec = [];
end

% True number of electrodes
nch = size(A,1);

% Number of permutations per size = n_perm

% What fractions of electrodes to take = e_f
e_n = nch-ceil(e_f*nch);
n_f = length(e_f); % number of electrodes to ignore

%% Ensure that all electrodes in loc are in adjacency matrix
% It should be perfectly aligned
for i = 1:length(pt(whichPt).new_elecs.electrodes)
    e_name = pt(whichPt).new_elecs.electrodes(i).name;
    found_it = 0;
    for j = 1:length(adj(7).data.labels)
        if strcmp(e_name,adj(7).data.labels{j}) == 1 && adj(7).data.ignore(j) == 0
            found_it = 1;
            break
        end
    end
    if found_it == 0
        error('Did not find %s in adjacency matrix.\n',e_name);
    end
end


%% Initialize matrices of resampled metrics

% Initialize cell array for each channel of the control centralities for
% each fraction and permutation
all_c_c = nan(nch,n_f,n_perm);
cc_reg = nan(nch,n_f,n_perm);
all_cc_norm = nan(nch,n_f,n_perm);

% Initialize cell array for each channel of the node strengths for each
% fraction and permutation
all_ns = nan(nch,n_f,n_perm);
all_ns_norm = nan(nch,n_f,n_perm);

% Initialize cell array for each channel of the betweenness centralities
% for each fraction and permutation
all_bc = nan(nch,n_f,n_perm);
all_bc_norm = nan(nch,n_f,n_perm);

% Initialize array for synhronizability for each fraction and permutation
all_sync = nan(n_f,n_perm);
all_sync_norm = nan(n_f,n_perm);

% Initialize array for efficiency for each fraction and permutation
all_eff = nan(n_f,n_perm);
all_eff_norm = nan(n_f,n_perm);

% Initialize array for eigenvector centrality
all_ec = nan(nch,n_f,n_perm);
all_ec_norm = nan(nch,n_f,n_perm);

% Initialize array for local efficiency
all_le = nan(nch,n_f,n_perm);

% Initialize array for clustering coefficient
all_clust = nan(nch,n_f,n_perm);
all_clust_norm = nan(nch,n_f,n_perm);

% Get arrays representing overlap and distance between removed channels and soz
overlap_soz = nan(n_f,n_perm);
dist_soz = nan(n_f,n_perm);
sz_soz_dist = nan(n_f,n_perm);

% Get arrays representing overlap and distance between removed channels and
% resected channels
overlap_resec = nan(n_f,n_perm);
dist_resec = nan(n_f,n_perm);
dist_nearest_resec = nan(n_f,n_perm);

% Initialize array representing electrodes in most synchronizig region
elecs_min = [];

% Participation coefficient
all_par = nan(nch,n_f,n_perm);

% Transitivity
all_trans = nan(n_f,n_perm);
all_trans_norm = nan(n_f,n_perm);

% Get true participation coefficient of electrodes (to determine the
% average participation coefficient of the removed electrodes)
[Ci,~]=modularity_und(A);
true_par = participation_coef(A,Ci);
avg_par_removed = nan(n_f,n_perm);

% Get true bc of electrodes
true_bc = betweenness_centrality(A,1);
avg_bc_removed = nan(n_f,n_perm);

%%  Loop through fractions
for f = 1:n_f
    
    
    %%  Do n_perm permutations per fraction
    for i_p = 1:n_perm
       
        if mod(i_p,100) == 0
            fprintf('Doing permutation %d of %d for removal fraction %d of %d...\n',...
                i_p,n_perm,f,n_f);
        end
        
        %% Take e_n electrodes and just remove them
        % (e_n is nch - n_f(f))
        if contig == 0 % random electrodes
            which_elecs = randperm(nch,e_n(f));
        elseif contig == 1 
            %{
            PRIMARY APPROACH
            % the e_n nearest neighbors, NOT RANDOM AT ALL because I am not
            % adding jitter and I am looping through channel by channel
            %}
            which_elecs = pickConChs(locs,e_n(f),0,0,i_p);
            
            
            % ALTERNATIVE APPROACH
            % Take just the individual channel
            % which_elecs = i_p;
        elseif contig == 2
            % Here I will take a random sample of all electrodes NOT
            % forming the SOZ, equal in number to number of SOZ elecs
            
            % Find non-soz electrodes
            n_soz = length(soz);
            soz_binary = ismember(1:nch,soz);
            not_soz = find(~soz_binary);
            
            % pick random sample equal in size to soz
            which_elecs = randsample(not_soz,min(n_soz,length(not_soz)));
            
        elseif contig == 3
            
            % just take soz
            which_elecs = soz;
        
        elseif contig == 4
            % Take random sample of ALL electrodes, equal in number to
            % number of soz elecs
            n_soz = length(soz);
            which_elecs = randperm(nch,n_soz);
            
        end
        
        %% Compare electrodes to SOZ and resection zone
        [overlap_soz_t,dist_soz_t,~] = compare_elecs(which_elecs,soz,locs,doPlot);
        [~,sz_soz_dist_t,~] = compare_elecs(which_elecs,sz_soz,locs,doPlot);
        % Get distance from SOZ and overlap with SOZ
        
        if isempty(resec) == 0
            [overlap_resec_t,dist_resec_t,dist_nearest_resec_t] = compare_elecs(which_elecs,resec,locs,doPlot);
            % Get distance from resection zone and overlap with resection
            % zone
        else
            overlap_resec_t = nan; dist_resec_t = nan; dist_nearest_resec_t = nan;
        end
        
        % Note that dist_resec(f,i_p) tells us the average distance between
        % the e_n nearest neighbors to channel i_p (because in the SOZ
        % analysis n_perm is the number of channels and we are looping
        % though each channel one by one) and the resection zone. We are
        % ignoring those e_n channels  
        overlap_soz(f,i_p) = overlap_soz_t;
        dist_soz(f,i_p) = dist_soz_t;
        overlap_resec(f,i_p) = overlap_resec_t;
        dist_resec(f,i_p) = dist_resec_t;
        dist_nearest_resec(f,i_p) = dist_nearest_resec_t;
        sz_soz_dist(f,i_p) = sz_soz_dist_t;
        
        %% Get bc and pc of removed electrodes
        avg_par_removed(f,i_p) = mean(true_par(which_elecs));
        avg_bc_removed(f,i_p) = mean(true_bc(which_elecs));
        
        %% Remove the electrodes from the new adjacency matrix
        % Note that this is one of the trickiest parts of this and the most
        % susceptible to error. I need to keep track of WHICH channels I am
        % reporting the updated nodal measures on.
        A_temp = A;
        ch_ids = 1:nch;
        A_temp(which_elecs,:) = [];
        A_temp(:,which_elecs) = [];
        ch_ids(which_elecs) = [];
        
        % Get new control centrality of remaining electrodes in the new
        % network
        c_c = control_centrality(A_temp);
        if i_p == 1
            cc_fake = nan(100,1);
            for i = 1:100
                cc_fake(i) = mean(control_centrality(generate_fake_graph(A_temp)));
            end
            m_cc_fake = mean(cc_fake);
        end
        
        
        % Get new node strength
        ns = node_strength(A_temp);
        if i_p == 1
            ns_fake = nan(100,1);
            for i = 1:100
                ns_fake(i) = mean(node_strength(generate_fake_graph(A_temp)));
            end
            m_ns_fake = mean(ns_fake);
        end
        
        % Get new clustering coefficient
        clust = clustering_coef_wu(A_temp);
        if i_p == 1
            clust_fake = nan(100,1);
            for i = 1:100
                clust_fake(i) = mean(clustering_coef_wu(generate_fake_graph(A_temp)));
            end
            m_clust_fake = mean(clust_fake);
        end
        
        % Get new betweenness centrality
        bc = betweenness_centrality(A_temp,1);
        if i_p == 1
            bc_fake = nan(100,1);
            for i = 1:100
                bc_fake(i) = mean(betweenness_centrality(generate_fake_graph(A_temp),1));
            end
            m_bc_fake = mean(bc_fake);
        end
        
        % Get new participation coefficient
        [Ci,~]=modularity_und(A_temp);
        par = participation_coef(A_temp,Ci);
        
        % get new synchronizability
        all_sync(f,i_p) = synchronizability(A_temp);
        
        % To normalize it, generate fake matrix just once, for the first
        % permutation, and use the same fake matrix each time. (Then
        % re-make it when drop to a new number of removed electrodes).
        if i_p == 1
            sync_fake = nan(100,1);
            for i = 1:100
                sync_fake(i) = synchronizability(generate_fake_graph(A_temp));
            end
            m_sync_fake = mean(sync_fake);
        end
        all_sync_norm(f,i_p) = synchronizability(A_temp)/m_sync_fake;
        
        % Get new efficiency
        all_eff(f,i_p) = efficiency_wei(A_temp, 0);
        if i_p == 1
            eff_fake = nan(100,1);
            for i = 1:100
                eff_fake(i) = efficiency_wei(generate_fake_graph(A_temp),0);
            end
            m_eff_fake = mean(eff_fake);
        end
        all_eff_norm(f,i_p) = efficiency_wei(A_temp,0)/m_eff_fake;
        
        % Get new local efficiency
        %le = efficiency_wei(A_temp,1);
        
        % get new transitivity
        all_trans(f,i_p) = transitivity_wu(A_temp);
        if i_p == 1
            trans_fake = nan(100,1);
            for i = 1:100
                trans_fake(i) = transitivity_wu(generate_fake_graph(A_temp));
            end
            m_trans_fake = mean(trans_fake);
        end
        all_trans_norm(f,i_p) = transitivity_wu(A_temp)/m_trans_fake;
        
        % get new eigenvector centrality
        ec = eigenvector_centrality_und(A_temp);
        if i_p == 1
            ec_fake = nan(100,1);
            for i = 1:100
                ec_fake(i) = mean(eigenvector_centrality_und(generate_fake_graph(A_temp)));
            end
            m_ec_fake = mean(ec_fake);
        end
        
        % new regional control centrality
        if isempty(pt(whichPt).resec) == 0 && e_f(f) > 0.2
            num_resec = length(pt(whichPt).resec.nums);
            if num_resec < length(ch_ids)
                [cc_regional,elecs_regional] = regional_control_centrality...
                    (A_temp,num_resec,locs(ch_ids,:),1);
                [~,min_cc_regional_true] = min(cc_regional);
                elecs_regional_min = elecs_regional(min_cc_regional_true,:);
                elecs_min = [elecs_min,ch_ids(elecs_regional_min)];
            else
                cc_regional = nan(length(ch_ids),1);
            end
            
            if 1 == 0
                % Plot the new min cc region and the electrodes we're
                % ignoring
                figure
                scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k')
                hold on
                scatter3(locs(ch_ids(elecs_regional_min),1),...
                    locs(ch_ids(elecs_regional_min),2),...
                    locs(ch_ids(elecs_regional_min),3),100,'b','filled')
                scatter3(locs(which_elecs,1),locs(which_elecs,2),locs(which_elecs,3),40,'k','filled')
                pause
                close(gcf)
            end
        else
            cc_regional = nan(length(ch_ids),1);
        end
        
        
        % Populate the node-specific metrics
        % If there were 10 electrodes to start, and I ignored 3 and 8, then
        % ch_ids would be [1 2 4 5 6 7 9 10]. When i is 3, then this
        % corresponds to ch_ids(3) = 4, and so I would be reporting nodal
        % metrics for original channel number 4 (which is channel number 3
        % in the resampled network).
        for i = 1:length(ch_ids)
            all_c_c(ch_ids(i),f,i_p) = c_c(i);
            all_ns(ch_ids(i),f,i_p) = ns(i);
            all_bc(ch_ids(i),f,i_p) = bc(i);
            all_par(ch_ids(i),f,i_p) = par(i);
            all_ec(ch_ids(i),f,i_p) = ec(i);
            all_clust(ch_ids(i),f,i_p) = clust(i);
            cc_reg(ch_ids(i),f,i_p) = cc_regional(i);
            %all_le(ch_ids(i),f,i_p) = le(i);
            
            all_cc_norm(ch_ids(i),f,i_p) = c_c(i)./m_cc_fake;
            all_ns_norm(ch_ids(i),f,i_p) = ns(i)./m_ns_fake;
            all_bc_norm(ch_ids(i),f,i_p) = bc(i)./m_bc_fake;
            all_ec_norm(ch_ids(i),f,i_p) = ec(i)./m_ec_fake;
            all_clust_norm(ch_ids(i),f,i_p) = clust(i)./m_clust_fake;
        end
        
        
    end
    
end

end