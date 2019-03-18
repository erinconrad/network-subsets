
%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);
outFolder = [resultsFolder,'../figures/'];

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);


% HUP078
whichPt = 8;
offset = [-3 27 7.7653]; % HUP078
%offset = [-1.0029,2.3087,28.9465]; %HUP068

figure
set(gcf,'Position',[175 369 500 500]);

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
p = plotGIFTI(g);
hold on
view(-120,-11);
alpha(p,0.4)
scatter3(new_locs(:,1),new_locs(:,2),new_locs(:,3),circSize,'k','filled');
plot3(new_locs([20,30],1),new_locs([20,30],2),new_locs([20,30],3),'k','linewidth',3)
print(gcf,[outFolder,'network_2'],'-depsc');

