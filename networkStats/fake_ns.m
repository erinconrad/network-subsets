function fake_ns(whichPt)

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);
outFolder = [resultsFolder,'fake_brains/'];

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);

% Load adj
[adj,~] = reconcileAdj(pt,whichPt,1);
A_all = adj(4).data;
A = squeeze(A_all(ceil(size(A_all,1)/2)-0,:,:));
ns = node_strength(A);

[~,all_ns] = resampleNetwork(A,10,[0.2 0.4 0.6 0.8 1],0,pt,whichPt,adj,1);
ns_rel = alt_rel_nodal(all_ns,ns);

locs = pt(whichPt).new_elecs.locs;
figure
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k')
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),100,ns,'filled')

for i = 1:size(all_ns,3)
figure
scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k')
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),100,squeeze(all_ns(:,1,i)),'filled')
pause
close(gcf)
end
end