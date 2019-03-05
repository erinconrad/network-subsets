function [all_c_c,all_ns,all_bc,all_sync,all_eff,overlap_soz,dist_soz,...
    overlap_resec,dist_resec,temp_centroid_min,elecs_min] = ...
    resampleNetwork(A,n_perm,e_f,contig,pt,whichPt,adj)

%{
This function resamples the network by removing a fraction of nodes and then
recalculates various network metrics    
%}
    
    
%% Get basic numbers

doPlot = 0; % Plot overlap of electrodes?

% Get electrode locs
locs = pt(whichPt).new_elecs.locs;

% Get soz electroes
soz = pt(whichPt).soz.nums;

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
n_f = length(e_f);

%% Remove from locs the electrodes that are not in the adjacency matrix
% These are electrodes that are in the location file but for which there is
% no EEG tracing

% THIS SHOULDN'T HAPPEN ANYMORE!
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

% Get arrays representing overlap and distance between removed channels and soz
overlap_soz = nan(n_f,n_perm);
dist_soz = nan(n_f,n_perm);

% Get arrays representing overlap and distance between removed channels and
% resected channels
overlap_resec = nan(n_f,n_perm);
dist_resec = nan(n_f,n_perm);

% Initialize array representing centroid of minimum
% control centrality region
temp_centroid_min =nan(n_f,n_perm,3);

elecs_min = [];

%%  Loop through fractions
for f = 1:n_f
    
    %%  Do n_perm permutations per fraction
    for i_p = 1:n_perm
        
        %% Take e_n electrodes and just remove them
        % (e_n is nch - n_f)
        if contig == 0 % random electrodes
            which_elecs = randperm(nch,e_n(f));
        elseif contig == 1 % random electrodes close to each other
            which_elecs = pickConChs(locs,e_n(f),0,20);
        end
        
        %% Compare electrodes to SOZ and resection zone
        [overlap_soz_t,dist_soz_t] = compare_elecs(which_elecs,soz,locs,doPlot);
        if isempty(resec) == 0
            [overlap_resec_t,dist_resec_t] = compare_elecs(which_elecs,resec,locs,doPlot);
        else
            overlap_resec_t = nan; dist_resec_t = nan;
        end
        overlap_soz(f,i_p) = overlap_soz_t;
        dist_soz(f,i_p) = dist_soz_t;
        overlap_resec(f,i_p) = overlap_resec_t;
        dist_resec(f,i_p) = dist_resec_t;
        
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
        all_sync(f,i_p) = synchronizability(A_temp)/...
            synchronizability(generate_fake_graph(A_temp));
        
        % Get new efficiency
        all_eff(f,i_p) = efficiency_wei(A_temp, 0)/...
            efficiency_wei(generate_fake_graph(A_temp),0);
        
        % new regional control centrality
        if isempty(pt(whichPt).resec) == 0 && e_f(f) == 0.8
            num_resec = length(pt(whichPt).resec.nums);
            [cc_regional,elecs_regional] = regional_control_centrality...
                (A_temp,num_resec,locs(ch_ids,:),1);
            [~,min_cc_regional_true] = min(cc_regional);
            elecs_regional_min = elecs_regional(min_cc_regional_true,:);
            elecs_min = [elecs_min,elecs_regional_min];
            temp_centroid_min(f,i_p,:) = mean(locs(ch_ids(elecs_regional_min),:));
            
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
            temp_centroid_min(f,i_p,:) = [nan nan nan];
        end
        
        
        % Populate the node-specific metrics
        for i = 1:length(ch_ids)
            all_c_c(ch_ids(i),f,i_p) = c_c(i);
            all_ns(ch_ids(i),f,i_p) = ns(i);
            all_bc(ch_ids(i),f,i_p) = bc(i);
        end
        
        
    end
    
end

end