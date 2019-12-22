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
[adj,~] = reconcileAdj(pt,whichPt,1);
A_all = adj(4).data;
A = squeeze(A_all(ceil(size(A_all,1)/2)-5,:,:));
cc = control_centrality(A);
all_cc = resampleNetwork(A,1,0.8,0,pt,whichPt,adj,1);
cc_diff = (all_cc-cc)./cc;

locs = pt(whichPt).new_elecs.locs;

% Load gifti
brainFolder = '/Users/erinconrad/Desktop/residency stuff/R25/actual work/data/brains/';
giftiFolder = [brainFolder,pt(whichPt).name,'/'];
names = dir([giftiFolder,'*pial.gii']);
fname2 = names(1).name;
g = gifti([giftiFolder,fname2]);

%% Plot network
if 1
figure
imagesc(A)
c = colorbar;
set(c,'ticklabels','')
set(c,'ytick',[])
set(gca,'xticklabels','')
set(gca,'yticklabels','')
set(gca,'xtick',[])
set(gca,'ytick',[])
print(gcf,[outFolder,'Fig1_network'],'-depsc');
close(gcf)
end

if 1
%% Plot resampled cc's
circSize = 250;
figure
%p = plotGIFTI(g);
%alpha(p,0.3)
hold on
view(-120,-11);
%alpha(p,0.4)
scatter3(new_locs(:,1),new_locs(:,2),new_locs(:,3),circSize+100,'k','filled');
hold on
scatter3(new_locs(~isnan(all_cc),1),new_locs(~isnan(all_cc),2),new_locs(~isnan(all_cc),3),...
    circSize,all_cc(~isnan(all_cc)),'filled');
scatter3(new_locs(isnan(all_cc),1),new_locs(isnan(all_cc),2),new_locs(isnan(all_cc),3),...
    circSize,'w','filled')
scatter3(new_locs(isnan(all_cc),1),new_locs(isnan(all_cc),2),new_locs(isnan(all_cc),3),...
    circSize,'rx','linewidth',3)
set(gca,'clim',prctile(cc,[10 90]));
set(gca,'fontsize',20)
xticklabels([])
yticklabels([])
zticklabels([])
grid off
axis off

% Get names of current figures
listing = dir(outFolder);
fname_start = 'resample_cc_no_brain';
all_nums = [];
for i = 1:length(listing)
    if contains(listing(i).name,fname_start)
        startIndex = regexp(listing(i).name,'\d+');
        endIndex = regexp(listing(i).name,'\.');
        num = str2double(listing(i).name(startIndex:endIndex-1));
        all_nums = [all_nums,num];
    end
end
if isempty(all_nums), curr_num = 1; else, curr_num = max(all_nums) + 1; end
print(gcf,[outFolder,fname_start,sprintf('%d',curr_num)],'-depsc');
close(gcf)

end

if 0
%% Plot removed electrodes
figure
p = plotGIFTI(g);
alpha(p,0.3)
hold on
view(-120,-11);
scatter3(new_locs(:,1),new_locs(:,2),new_locs(:,3),circSize+100,'k','filled');
scatter3(new_locs(:,1),new_locs(:,2),new_locs(:,3),...
    circSize,'w','filled');
scatter3(new_locs(isnan(all_cc),1),new_locs(isnan(all_cc),2),new_locs(isnan(all_cc),3),...
    circSize,'rx','linewidth',3)
print(gcf,[outFolder,'Fig1_rm'],'-depsc');
close(gcf)

%% Plot original cc
figure
p = plotGIFTI(g);
alpha(p,0.3)
hold on
view(-120,-11);
scatter3(new_locs(:,1),new_locs(:,2),new_locs(:,3),circSize+100,'k','filled');
scatter3(new_locs(:,1),new_locs(:,2),new_locs(:,3),...
    circSize,cc,'filled');
set(gca,'clim',prctile(cc,[10 90]));
print(gcf,[outFolder,'Fig1_original'],'-depsc');
close(gcf)
end

end
