function save_ns(whichPts)

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);

% High gamma, EEC
freq = 'high_gamma';
which_sec = 0;


out_folder = [resultsFolder,'resec_ns/'];

if exist(out_folder,'dir') == 0
    mkdir(out_folder);
end


for whichPt = whichPts
    
    name = pt(whichPt).name;
    
    % Skip if all electrode locations are -1 (means we don't have electrode
    % locations)
    if unique(pt(whichPt).new_elecs.locs) == -1
        continue
    end
    
    % Skip if no resection
    if isempty(pt(whichPt).resec) == 0
        resec = pt(whichPt).resec.nums;
    else
        continue
    end
        
    ch_nums = 1:length(pt(whichPt).new_elecs.electrodes);
    yes_resec = ismember(ch_nums,resec);
    no_resec = ~yes_resec;
    
    %% Get adjacency matrix
    [adj,~] = reconcileAdj(pt,whichPt,1);
    if isempty(adj) == 1
        fprintf('Cannot do %s\n\n',name);
        continue;
    end

    %% Get appropriate frequency band
    if strcmp(freq,'high_gamma') == 1
        A_all = adj(4).data;
        if contains(adj(4).name,'highgamma') == 0
            error('This isn''t gamma!'\n');
        end
    end
    
    
    %% Get which time
    % Start with the middle and add which second
    if ceil(size(A_all,1)/2)+which_sec <= 0, continue; end
    if ceil(size(A_all,1)/2)+which_sec > size(A_all,1), continue; end
    A = squeeze(A_all(ceil(size(A_all,1)/2)+which_sec,:,:));
    if sum(sum(isnan(A))) == sum(sum(ones(size(A))))
        continue
    end
    
    %% Get node strength
    ns = node_strength(A);
    
    %% Save ns into a struct to do more stats
    ns(whichPt).name = name;
    ns(whichPt).ns = ns;
    ns(whichPt).yes_resec = yes_resec;
    ns(whichPt).no_resec = no_resec;
    
    
end

%% Save the structure
save([out_folder,'ns.mat'],'ns');