function methods_fig

% I will look at control centrality for a single patient. I will then
% randomly remove 20% of the nodes and recalculate control centrality. I
% will then show the correlation for each electrode of old and new control
% centrality (3 subplots)

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);
outFolder = [resultsFolder,'figures/'];

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

figure
set(gcf,'Position',[175 369 1300 910]);
[ha,pos] = tight_subplot(2,2,[0.05 0.03],[0.01 0.01],[0.02 0.03]);


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

%% Plot original cc's
axes(ha(1));
p = plotGIFTI(g);
hold on
view(-120,-11);
%alpha(p,0.4)
scatter3(new_locs(:,1),new_locs(:,2),new_locs(:,3),circSize,'k');
scatter3(new_locs(:,1),new_locs(:,2),new_locs(:,3),circSize,cc,'filled');
set(gca,'clim',prctile(cc,[10 90]));
colorbar
annotation('textbox',[0.10 0.9 0.1 0.1],'String',...
    'Electrode control centrality','LineStyle','none','fontsize',25);
set(gca,'fontsize',20)

%% Plot original cc's but indicate which electrodes are removed
axes(ha(2));
p = plotGIFTI(g);
hold on
view(-120,-11);
%alpha(p,0.4)
scatter3(new_locs(:,1),new_locs(:,2),new_locs(:,3),circSize,'k');
scatter3(new_locs(:,1),new_locs(:,2),new_locs(:,3),circSize,'w','filled');
%{
scatter3(new_locs(isnan(all_cc),1),new_locs(isnan(all_cc),2),new_locs(isnan(all_cc),3),...
    circSize,'r','filled');
%}
scatter3(new_locs(isnan(all_cc),1),new_locs(isnan(all_cc),2),new_locs(isnan(all_cc),3),...
    circSize+20,'x','r','linewidth',2);
annotation('textbox',[0.64 0.9 0.1 0.1],'String',...
    'Ignored electrodes','LineStyle','none','fontsize',25);
set(gca,'fontsize',20)


%% Plot resampled cc's
axes(ha(3));
p = plotGIFTI(g);
hold on
view(-120,-11);
%alpha(p,0.4)
scatter3(new_locs(~isnan(all_cc),1),new_locs(~isnan(all_cc),2),new_locs(~isnan(all_cc),3),circSize,'k');
scatter3(new_locs(~isnan(all_cc),1),new_locs(~isnan(all_cc),2),new_locs(~isnan(all_cc),3),...
    circSize,all_cc(~isnan(all_cc)),'filled');
set(gca,'clim',prctile(cc,[10 90]));
colorbar
annotation('textbox',[0.11 0.4 0.1 0.1],'String',...
    'New control centrality','LineStyle','none','fontsize',25);
set(gca,'fontsize',20)

%% Plot difference
axes(ha(4));
p = plotGIFTI(g);
hold on
view(-120,-11);
%alpha(p,0.4)
scatter3(new_locs(~isnan(all_cc),1),new_locs(~isnan(all_cc),2),new_locs(~isnan(all_cc),3),circSize,'k');
scatter3(new_locs(~isnan(all_cc),1),new_locs(~isnan(all_cc),2),new_locs(~isnan(all_cc),3),...
    circSize,cc_diff(~isnan(all_cc)),'filled');
set(gca,'clim',prctile(cc_diff,[10 90]));
colorbar
annotation('textbox',[0.63 0.4 0.1 0.1],'String',...
    'Relative change','LineStyle','none','fontsize',25);
set(gca,'fontsize',20)

saveas(gcf,[outFolder,'Fig1.fig'])
print(gcf,[outFolder,'Fig1'],'-depsc');
saveas(gcf,[outFolder,'Fig1.eps'])
saveas(gcf,[outFolder,'Fig1.png'])


end