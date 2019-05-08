function test_nodal_size

sizes = [20 40 60 80];
np = 100;

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

ns_rel_old = zeros(length(sizes),1);
cc_rel_old = zeros(length(sizes),1);
ns_rel_new = zeros(length(sizes),1);
cc_rel_new = zeros(length(sizes),1);


% Get original
ns_orig = node_strength(fake_A_certain_size(100));
cc_orig = control_centrality(fake_A_certain_size(100));


for s = 1:length(sizes)
    
    fprintf('Doing size %d of %d\n',s,length(sizes))
    
    count = 0;
    
    ns_all = zeros(np,sizes(s));
    cc_all = zeros(np,sizes(s));
    
    for p = 1:np
      
        
        % Make a fake graph
        A = fake_A_certain_size(sizes(s));
        
        % Calculate mean of measures
        ns = (node_strength(A));
        c_c = (control_centrality(A));
        bc = (betweenness_centrality(A,1));
        ec = (eigenvector_centrality_und(A));
        clust = (clustering_coef_wu(A));
        
       
        ns_all(p,:) = ns;
        cc_all(p,:) = c_c;
        
    
    end
    

    % Old way of reliability
    ns_rel_old(s) = nanstd(ns_orig).^2/(nanstd(ns_orig).^2+mean(nanstd(ns_all,0,1).^2));
    cc_rel_old(s) = nanstd(cc_orig).^2/(nanstd(cc_orig).^2+mean(nanstd(cc_all,0,1).^2));
    
    % new way of reliability
    ns_rel_new(s) = nanmean(nanstd(ns_all,0,2).^2)/(nanmean(nanstd(ns_all,0,2).^2)+mean(nanstd(ns_all,0,1).^2));
    cc_rel_new(s) = nanmean(nanstd(cc_all,0,2).^2)/(nanmean(nanstd(cc_all,0,2).^2)+mean(nanstd(cc_all,0,1).^2));
    
end


end