function aes_network

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
addpath(genpath(scriptFolder));
outFolder = [resultsFolder,'figures/'];

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);


% HUP078
whichPt = 8;

% Load adj
[adj,~] = reconcileAdj(pt,whichPt,1);
A_all = adj(4).data;
A = squeeze(A_all(ceil(size(A_all,1)/2)-5,:,:));


figure
set(gcf,'Position',[175 369 910 910]);
tight_subplot(1, 1, [0.1 0.04], [0.13 0.08],[0.14 0.03]);
imagesc(A)
xlabel('Electrode number')
ylabel('Electrode number')
title('Connection strength')
%{
hcb = colorbar;
title(hcb,'Connection strength')
set(hcb,'ticklabels',[])
%}
set(gca,'fontsize',40)

end