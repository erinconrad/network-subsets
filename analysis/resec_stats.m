function resec_stats

% This just compares average n_s for resection zone to other areas
doPlot = 0;

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

%% Load ns structure
load([resultsFolder,'resec_ns/ns.mat']);
load([dataFolder,'structs/info.mat']);

all_resec_ns = [];
all_no_resec_ns = [];

for i = 1:length(ns)
    
    if isempty(ns(i).ns) == 1, continue; end
    
    ns_temp = ns(i).ns;
    yes_resec = ns(i).yes_resec;
    no_resec = ns(i).no_resec;
    
    %% Plot
    if doPlot == 1
        locs = pt(i).new_elecs.locs;
        scatter3(locs(:,1),locs(:,2),locs(:,3),200,'k');
        hold on
        scatter3(locs(:,1),locs(:,2),locs(:,3),200,ns_temp,'filled');
        scatter3(locs(yes_resec,1),locs(yes_resec,2),locs(yes_resec,3),70,'k','filled');
        colorbar
        pause
        close(gcf)
    end
    
    %% Get average ns of resection zone versus outside it
    ns_resec = mean(ns_temp(yes_resec));
    ns_no_resec = mean(ns_temp(no_resec));
    
    all_resec_ns = [all_resec_ns;ns_resec];
    all_no_resec_ns = [all_no_resec_ns;ns_no_resec];
    
    
end

%% Descriptive stats
fprintf(['The mean node strength of electrodes in the resection zone is\n'...
    '%1.2f (range %1.2f-%1.2f), and the mean node strength of electrodes\n'...
    'not in the resection zone is %1.2f (range %1.2f-%1.2f)\n'],mean(all_resec_ns),...
    min(all_resec_ns), max(all_resec_ns), mean(all_no_resec_ns),...
    min(all_no_resec_ns), max(all_no_resec_ns));

%% Two sided paired t test
[~,p,~,stats] = ttest(all_resec_ns,all_no_resec_ns);
fprintf(['Paired t-test: t = %1.3f, p = %1.3f\n'],stats.tstat,p)


end