
%% Parameters
np = 10;
e_f = 0.8; 
freq = 'high_gamma';
contigs = 0;
which_sec = -5;
whichPt = 8;
circSize = 100;
delay = 1; % time delay between steps

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);

name = pt(whichPt).name;
fsave = [resultsFolder,'../','figures/min_cc_movie/movie1.gif'];
fsave_image = [resultsFolder,'../','figures/min_cc_movie/im1'];

%% Get adjacency matrix
[adj,~] = reconcileAdj(pt,whichPt);

% Get electrode locs
locs = pt(whichPt).new_elecs.locs;
nchs = size(locs,1);
e_n = nchs-ceil(e_f*nchs);
num_resec = length(pt(whichPt).resec.nums);

if strcmp(freq,'high_gamma') == 1
    A_all = adj(4).data;
    if contains(adj(4).name,'highgamma') == 0
        error('This isn''t gamma!'\n');
    end
end

% Start with the middle and add which second
A = squeeze(A_all(ceil(size(A_all,1)/2)+which_sec,:,:));
if nchs ~= size(A,1), error('what\n'); end

%% Get true regional cc
% Get regional control centralities
[cc_regional1,elecs_regional] = regional_control_centrality(A,num_resec,locs,1);

% Get identity of region with lowest regional control centrality
[~,min_cc_regional_true] = min(cc_regional1);

% Get electrodes in region with lowest cc
elecs_regional_min_true = elecs_regional(min_cc_regional_true,:);

%% Loop through permutations and recalculate regional cc
min_elecs_perm = cell(np+1,1);
min_elecs_perm{1} = elecs_regional_min_true;
all_elecs_perm = cell(np+1,1);
all_elecs_perm{1} = 1:nchs;
for ip = 2:np+1
    
    % Select random group of channels and remove them
    which_elecs = randperm(nchs,e_n);
    A_temp = A;
    ch_ids = 1:nchs;
    A_temp(which_elecs,:) = [];
    A_temp(:,which_elecs) = [];
    ch_ids(which_elecs) = [];
    all_elecs_perm{ip} = ch_ids;
    
    % Calculate new regional cc
    [cc_regional,elecs_regional] = regional_control_centrality...
        (A_temp,num_resec,locs(ch_ids,:),1);
    [~,min_cc_regional_true] = min(cc_regional);
    elecs_regional_min = elecs_regional(min_cc_regional_true,:);
    min_elecs_perm{ip} = ch_ids(elecs_regional_min);
end

%% Image
figure
set(gcf,'position',[440 405 888 393]);
[ha,~] = tight_subplot(1,2,[0.01 0.03],[0.01 0.09],[0.01 0.01]);
axes(ha(1))
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k');
hold on
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,cc_regional1,'filled');
xticklabels([]);
yticklabels([]);
zticklabels([]);
title('Regional control centrality');
set(gca,'fontsize',20)
    
axes(ha(2))
scatter3(locs(:,1),locs(:,2),locs(:,3),circSize,'k');
hold on
min_chs = min_elecs_perm{1};
scatter3(locs(min_chs,1),locs(min_chs,2),locs(min_chs,3),circSize,'r','filled');
xticklabels([]);
yticklabels([]);
zticklabels([]);
title('Most synchronizing region');
set(gca,'fontsize',20)
print(gcf,fsave_image,'-depsc');


%% Make movie
for i = 1:np
    fig = figure;
    
    chs = all_elecs_perm{i};
    min_chs = min_elecs_perm{i};
    scatter3(locs(chs,1),locs(chs,2),locs(chs,3),circSize,'k');
    hold on
    scatter3(locs(min_chs,1),locs(min_chs,2),locs(min_chs,3),circSize,'r','filled');
    xticklabels([]);
    yticklabels([]);
    zticklabels([]);
    
    if i == 1
        title('Most synchronizing region');
    else
        title(sprintf('Permutation %d',i-1));
    end
    set(gca,'fontsize',20);
    
     % capture the figure as a frame in the gif
    F(i) = getframe(fig);
    im = frame2im(F(i));
    [imind,cm] = rgb2ind(im,256);
    
    if i == 1
        imwrite(imind,cm,fsave,'gif', 'Loopcount',inf,'DelayTime',delay);
    else
        imwrite(imind,cm,fsave,'gif','WriteMode','append','DelayTime',delay);
    end

    close(fig)
    
end

