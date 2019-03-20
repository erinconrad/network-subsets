function compare_cc(true_cc_all,whichPt)

which_sec = 0;

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);

[adj,~] = reconcileAdj(pt,whichPt);
A_all = adj(4).data;
A = squeeze(A_all(ceil(size(A_all,1)/2)+which_sec,:,:));

true_cc = true_cc_all(ceil(size(true_cc_all,1)/2)+which_sec,:);

cc = control_centrality(A);

locs = pt(whichPt).new_elecs.locs; 

figure
subplot(1,2,1)
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
scatter3(locs(:,1),locs(:,2),locs(:,3),100,cc,'filled');
set(gca,'clim',prctile(cc,[10 90]));

subplot(1,2,2)
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
scatter3(locs(:,1),locs(:,2),locs(:,3),100,true_cc,'filled');
set(gca,'clim',prctile(true_cc,[10 90]));


end