function fake_cc

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);
outFolder = [resultsFolder,'fake_brains/'];

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);

% HUP078
whichPt = 8;
offset = [-3 27 7.7653]; % HUP078
%offset = [-1.0029,2.3087,28.9465]; %HUP068

% Load adj
[adj,~] = reconcileAdj(pt,whichPt);
A_all = adj(4).data;
A = squeeze(A_all(ceil(size(A_all,1)/2)-5,:,:));
cc = control_centrality(A);
all_cc = resampleNetwork(A,1,0.8,0,pt,whichPt,adj);
cc_diff = (all_cc-cc)./cc;

locs = pt(whichPt).new_elecs.locs;
A = transform_elecs_to_brain(pt,whichPt);
new_locs = A*locs-offset;

% Load gifti
brainFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/brains/';
giftiFolder = [brainFolder,pt(whichPt).name,'/'];
names = dir([giftiFolder,'*pial.gii']);
fname2 = names(1).name;
g = gifti([giftiFolder,fname2]);

circSize = 200;

figure
%% Plot resampled cc's
p = plotGIFTI(g);
hold on
view(-120,-11);
%alpha(p,0.4)
scatter3(new_locs(~isnan(all_cc),1),new_locs(~isnan(all_cc),2),new_locs(~isnan(all_cc),3),circSize,'k');
scatter3(new_locs(~isnan(all_cc),1),new_locs(~isnan(all_cc),2),new_locs(~isnan(all_cc),3),...
    circSize,all_cc(~isnan(all_cc)),'filled');
set(gca,'clim',prctile(cc,[10 90]));
colorbar
%{
annotation('textbox',[0.11 0.4 0.1 0.1],'String',...
    'New control centrality','LineStyle','none','fontsize',25);
%}
set(gca,'fontsize',20)
print(gcf,[outFolder,'Fig1_temp'],'-depsc');
close(gcf)

end
