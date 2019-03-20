function compare_cc(whichPts)

%{
The purpose of this is to compare the cc's I calculated from the ones John
calculated. When I ran this on ~25 patients they were the same.
%}

which_sec = -5;

if isempty(whichPts) == 1
    whichPts = [1:5,8:33];
end

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);
true_cc_folder_base = [dataFolder,'control_centralities/'];

for whichPt = whichPts

    % Load new true cc    
    name = pt(whichPt).name;
    true_cc_folder = [true_cc_folder_base,name,'/'];
    load([true_cc_folder,'cc_true.mat']);
    true_cc_all = cc;

    [adj,~] = reconcileAdj(pt,whichPt);
    A_all = adj(4).data;
    A = squeeze(A_all(ceil(size(A_all,1)/2)+which_sec,:,:));

    true_cc = true_cc_all(ceil(size(true_cc_all,1)/2)+which_sec,:)';
    %A = squeeze(A_all(1,:,:));
    %true_cc = true_cc_all(1,:)';

    cc_new = control_centrality(A);

    locs = pt(whichPt).new_elecs.locs; 
    
    max_dev = max(abs(cc_new-true_cc));
    if max_dev > 0.0001
        fprintf('WARNING, discrepancy for %s\n\n\n\n',name);
    else
        fprintf('Max dev for %s is %1.2e\n\n',name,max_dev);
    end
    
    
    figure
    subplot(1,2,1)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,cc_new,'filled');
    set(gca,'clim',prctile(cc_new,[10 90]));
    colorbar

    subplot(1,2,2)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,true_cc,'filled');
    set(gca,'clim',prctile(true_cc,[10 90]));
    colorbar
    %}


end




end