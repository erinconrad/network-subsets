function basic_info_networks(pt,stats)

%% Get number of patients with electrode info
np = 0;
for i = 1:length(stats)
    if isempty(stats(i).name) == 0
        np = np + 1;
    end
end

fprintf('There are %d patients with electrode info.\n\n',np);

%% Get number of patients with >1 seizure
np_sz_2 = 0;
for i = 1:length(stats)
    if isempty(stats(i).name) == 1, continue; end
    if strcmp(stats(i).name,pt(i).name) == 0, error('what\n'); end
    if length(pt(i).sz) == 1
        fprintf('%s has only one seizure\n',stats(i).name);
    elseif isempty(pt(i).sz) == 1
        error('what\n');
    else
        np_sz_2 = np_sz_2 + 1;
    end
end
fprintf('There are %d patients with >1 sz.\n\n',np_sz_2);

%% Get number of patients with resection info
np_resec = 0;
for i = 1:length(stats)
    if isempty(stats(i).name) == 1, continue; end
    if strcmp(stats(i).name,pt(i).name) == 0, error('what\n'); end
    
    if isempty(pt(i).resec) == 1
        fprintf('%s has no resection info\n',stats(i).name);
    else
        np_resec = np_resec + 1;
    end
    
end

fprintf('There are %d patients with resection info.\n\n',np_resec);

%% Get electrode info
n_elecs = [];
for i = 1:length(stats)
    n_elecs = [n_elecs;length(pt(i).new_elecs.electrodes)];
end

    

end