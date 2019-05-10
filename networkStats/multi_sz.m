function multi_sz(whichPts)

%% Parameters
e_f = [0.2 0.4 0.6 0.8 1];
nperm = 1e3;
freq = 4; %high-gamma
which_sec = 0; % EEC

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);

if exist([resultsFolder,'sz_comp.mat'],'file') ~= 0
    stats = load([resultsFolder,'sz_comp.mat'],'file');
    stats = stats.stats;
else
    stats = struct;
end

for whichPt = whichPts
    
    if unique(pt(whichPt).new_elecs.locs) == -1
        continue
    end
   
    stats(whichPt).name = pt(whichPt).name;
    
    % Skip if already done
    if length(stats) >= whichPt
        if isfield(stats(whichPt),'bc') == 1
            if isfield(stats(whichPt).bc,'sz') == 1
                if isfield(stats(whichPt).bc.sz(2),'perm') == 1
                    fprintf('Did %s, skipping\n',pt(whichPt).name);
                    continue
                end
            end
        end
    end
    fprintf('Doing %s removal percentage ',pt(whichPt).name);
    
    %% Get adjacency matrix
    adj1 = reconcileAdj(pt,whichPt,1);
    adj2 = reconcileAdj(pt,whichPt,2);
    if isempty(adj1) == 1 || isempty(adj2) == 1
        fprintf('Cannot do %s\n\n',name);
        continue;
    end
    
    %% Get the high gamma frequency
    A1 = adj1(freq).data;
    A2 = adj2(freq).data;
    
    %% Get EEC
    A1 = squeeze(A1(ceil(size(A1,1)/2)+which_sec,:,:));
    A2 = squeeze(A2(ceil(size(A2,1)/2)+which_sec,:,:));
    A(1,:,:) = A1;
    A(2,:,:) = A2;
    
    % True number of electrodes
    nchs = size(A1,1);
    
    % What fractions of electrodes to take = e_f
    e_n = nchs-ceil(e_f*nchs);
    n_f = length(e_f); % number of removal percentages
    
    %% Get metrics
    for s = 1:2
        A_temp = squeeze(A(s,:,:));
        
        % control centrality
        stats(whichPt).cc.sz(s).true = control_centrality(A_temp);
        
        % node strength
        stats(whichPt).ns.sz(s).true = node_strength(A_temp);
        
        % eigenvector centrality
        stats(whichPt).ec.sz(s).true = eigenvector_centrality_und(A_temp);
        
        % clustering coefficient
        stats(whichPt).clust.sz(s).true = clustering_coef_wu(A_temp);
        
        % betweenness centrality
        stats(whichPt).bc.sz(s).true = betweenness_centrality(A_temp,1);
        
        % synchronizability
        stats(whichPt).sync.sz(s).true = synchronizability(A_temp);
        
        % global efficiency
        stats(whichPt).eff.sz(s).true = efficiency_wei(A_temp,0);
        
        % transitivity
        stats(whichPt).trans.sz(s).true = transitivity_wu(A_temp);
        
    end
    
    %% Initialize permutation arrays
    for s = 1:2
        stats(whichPt).cc.sz(s).perm = nan(n_f,nperm,nchs);
        stats(whichPt).ns.sz(s).perm = nan(n_f,nperm,nchs);
        stats(whichPt).ec.sz(s).perm = nan(n_f,nperm,nchs);
        stats(whichPt).clust.sz(s).perm = nan(n_f,nperm,nchs);
        stats(whichPt).bc.sz(s).perm = nan(n_f,nperm,nchs);
        
        stats(whichPt).sync.sz(s).perm = nan(n_f,nperm);
        stats(whichPt).sync.sz(s).perm = nan(n_f,nperm);
        stats(whichPt).sync.sz(s).perm = nan(n_f,nperm);
        stats(whichPt).sync.sz(s).perm = nan(n_f,nperm);
    end
    
    %% Resample network
    for f = 1:n_f
        
        fprintf(' %d',f);
        
        for i_p = 1:nperm
        
            % Select electrodes for removal
            which_elecs = randperm(nchs,e_n(f));

            % Remove electrodes
            A_perm = A;
            ch_ids = 1:nchs;
            A_perm(:,which_elecs,:) = [];
            A_perm(:,:,which_elecs) = [];
            ch_ids(which_elecs) = [];

            % Loop through seizures
            for s = 1:2

                A_temp = squeeze(A_perm(s,:,:));
                
                cc = control_centrality(A_temp);
                ns = node_strength(A_temp);
                ec = eigenvector_centrality_und(A_temp);
                clust = clustering_coef_wu(A_temp);
                bc = betweenness_centrality(A_temp,1);
                sync = synchronizability(A_temp);
                eff = efficiency_wei(A_temp,0);
                trans = transitivity_wu(A_temp);
                
                
                % Get global metrics
                % synchronizability
                stats(whichPt).sync.sz(s).perm(f,i_p) = sync;

                % global efficiency
                stats(whichPt).eff.sz(s).perm(f,i_p) = eff;

                % transitivity
                stats(whichPt).trans.sz(s).perm(f,i_p) = trans;
                
                % Get nodal metrics
                for i = 1:length(ch_ids)

                    % control centrality
                    stats(whichPt).cc.sz(s).perm(f,i_p,ch_ids(i)) = cc(i);
                    
                    % node strength
                    stats(whichPt).ns.sz(s).perm(f,i_p,ch_ids(i)) = ns(i);
                    
                    % eigenvector centrality
                    stats(whichPt).ec.sz(s).perm(f,i_p,ch_ids(i)) = ec(i);
                    
                    % clustering coefficient
                    stats(whichPt).clust.sz(s).perm(f,i_p,ch_ids(i)) = clust(i);
                    
                    % betweenness centrality
                    stats(whichPt).bc.sz(s).perm(f,i_p,ch_ids(i)) = bc(i);                       
                    
                end

            end
        
        
        end
        
    end
    fprintf('\n');
    save([resultsFolder,'sz_comp.mat'],'stats');
    
end


end