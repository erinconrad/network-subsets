function all_c_c = resampleNetwork(A,n_perm,e_f,contig,locs)

% True number of electrodes
nch = size(A,1);

% Number of permutations per size = n_perm

% What fractions of electrodes to take = e_f
e_n = nch-ceil(e_f*nch);
n_f = length(e_f);

% Initialize cell array for each channel of the control centralities for
% each fraction and permutation
all_c_c = zeros(nch,n_f,n_perm);

% Initialize cell array for each channel of the node strengths for each
% fraction and permutation
all_ns = zeros(nch,n_f,n_perm);

% Initialize cell array for each channel of the betweenness centralities
% for each fraction and permutation
all_bc = zeros(nch,n_f,n_perm);

% Loop through fractions
for f = 1:n_f
    
    % Do n_perm permutations per fraction
    for i_p = 1:n_perm
        
        % Take e_n electrodes and just remove them
        
        if contig == 0
            which_elecs = randperm(nch,e_n(f));
        elseif contig == 1
            which_elecs = pickConChs(locs,e_n(f),0);
        end
        
        A_temp = A;
        ch_ids = 1:nch;
        A_temp(which_elecs,:) = [];
        A_temp(:,which_elecs) = [];
        ch_ids(which_elecs) = [];
        
        % Get control centrality
        c_c = control_centrality(A_temp);
        
        % Get node strength
        ns = node_strength(A_temp);
        
        % Get betweenness centrality
        bc = betweenness_centrality(A_temp);
        
        % MAKE SURE THAT I INDEXED THE CH ID CORRECTLY
        for i = 1:length(ch_ids)
            all_c_c(ch_ids(i),f,i_p) = c_c(i);
            all_ns(ch_ids(i),f,i_p) = ns(i);
            all_bc(ch_ids(i),f,i_p) = bc(i);
        end
        
        % fill in removed channels with nans
        for i = 1:length(which_elecs)
            all_c_c(which_elecs(i),f,i_p) = nan;
        end
        
    end
    
end

end