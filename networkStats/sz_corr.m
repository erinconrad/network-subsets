function sz_corr(whichPts)

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

for whichPt = whichPts
    
    if unique(pt(whichPt).new_elecs.locs) == -1
        continue
    end
    
    
     %% Get adjacency matrix
    adj1 = reconcileAdj(pt,whichPt,1);
    adj2 = reconcileAdj(pt,whichPt,2);
    if isempty(adj1) == 1 || isempty(adj2) == 1
        fprintf('Cannot do %s\n\n',pt(whichPt).name);
        continue;
    end
    
    %% Get the high gamma frequency
    A1 = adj1(freq).data;
    A2 = adj2(freq).data;
    
    %% Get EEC
    A = [];
    A1 = squeeze(A1(ceil(size(A1,1)/2)+which_sec,:,:));
    A2 = squeeze(A2(ceil(size(A2,1)/2)+which_sec,:,:));
    A(1,:,:) = A1;
    A(2,:,:) = A2;
    
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
       
        
    end
    
    fprintf('For %s, correlation between sz 1 and sz 2 is...\n',pt(whichPt).name);
    
    metric_names = {'cc','ns','ec','clust','bc'};
    for i =1:length(metric_names)
        base = stats(whichPt).(metric_names{i});
        
        rho = corr(base.sz(1).true,base.sz(2).true,'Type','Spearman');
        
        fprintf('%1.2f for %s,\n',rho,metric_names{i});
        
    end
    
    fprintf('\n\n');
end

end