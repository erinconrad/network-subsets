function [all_c_c,all_ns,all_bc,all_sync,all_eff] = ...
    resampleNetwork(A,n_perm,e_f,contig,locs)

%{
This function resamples the network by taking a fraction of nodes and then
recalculates various network metrics    
%}
    
    
%% Get basic numbers

% True number of electrodes
nch = size(A,1);

% Number of permutations per size = n_perm

% What fractions of electrodes to take = e_f
e_n = nch-ceil(e_f*nch);
n_f = length(e_f);

%% Initialize matrices of resampled metrics

% Initialize cell array for each channel of the control centralities for
% each fraction and permutation
all_c_c = nan(nch,n_f,n_perm);

% Initialize cell array for each channel of the node strengths for each
% fraction and permutation
all_ns = nan(nch,n_f,n_perm);

% Initialize cell array for each channel of the betweenness centralities
% for each fraction and permutation
all_bc = nan(nch,n_f,n_perm);

% Initialize array for synhronizability for each fraction and permutation
all_sync = nan(n_f,n_perm);

% Initialize array for efficiency for each fraction and permutation
all_eff = nan(n_f,n_perm);

%%  Loop through fractions
for f = 1:n_f
    
    %%  Do n_perm permutations per fraction
    for i_p = 1:n_perm
        
        %% Take e_n electrodes and just remove them
        % (e_n is nch - n_f)
        if contig == 0 % random electrodes
            which_elecs = randperm(nch,e_n(f));
        elseif contig == 1 % random electrodes close to each other
            which_elecs = pickConChs(locs,e_n(f),0);
        end
        
        % Remove the electrodes from the new adjacency matrix
        A_temp = A;
        ch_ids = 1:nch;
        A_temp(which_elecs,:) = [];
        A_temp(:,which_elecs) = [];
        ch_ids(which_elecs) = [];
        
        % Get new control centrality
        c_c = control_centrality(A_temp);
        
        % Get new node strength
        ns = node_strength(A_temp);
        
        % Get new betweenness centrality
        bc = betweenness_centrality(A_temp,1);
        
        % get new synchronizability
        all_sync(f,i_p) = synchronizability(A_temp);
        
        % Get new efficiency
        all_eff(f,i_p) = efficiency_wei(A_temp, 0);
        
        % Populate the node-specific metrics
        for i = 1:length(ch_ids)
            all_c_c(ch_ids(i),f,i_p) = c_c(i);
            all_ns(ch_ids(i),f,i_p) = ns(i);
            all_bc(ch_ids(i),f,i_p) = bc(i);
        end
        
        
    end
    
end

end