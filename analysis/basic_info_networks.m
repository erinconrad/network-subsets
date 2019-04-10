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

%% Get gender info
n_men = 0;
n_women = 0;
for i = 1:length(stats)
    if isempty(stats(i).name) == 1, continue; end
    if strcmp(pt(i).name,'HUP111B') == 1
        continue
    end
    
    if strcmp(pt(i).clinical.sex,'M') == 1
        n_men = n_men + 1;
    elseif strcmp(pt(i).clinical.sex,'F') == 1
        n_women = n_women + 1;
    else
        error('what\n');
    end
end

fprintf('There were %d women and %d men.\n',n_women,n_men);

%% Get age at surgery info
age = [];
for i = 1:length(stats)
    if isempty(stats(i).name) == 1, continue; end
    if strcmp(pt(i).name,'HUP111B') == 1
        continue
    end
    
    t = pt(i).clinical.ageSurgery;
    
    if contains(t,'-') == 1
        a = regexp(t,'-');
        num1 = str2num(t(1:a-1));
        num2 = str2num(t(a+1:end));
        age_num = mean([num1,num2]);
    elseif contains(t,'?') == 1
        age_num = nan;
    else
        age_num = str2num(t);
    end
    age = [age;age_num];
end

fprintf('The average age was %1.1f (range %1.1f-%1.1f).\n',nanmean(age),...
    min(age),max(age));

%% Get electrode info
n_elecs = [];
for i = 1:length(stats)
    if isempty(stats(i).name) == 1, continue; end
    n_elecs = [n_elecs;length(pt(i).new_elecs.electrodes)];
end

fprintf('The mean number of electrodes was %1.1f (range %d-%d across patients).\n',...
    mean(n_elecs),min(n_elecs),max(n_elecs));

%% Get n seizures
n_sz = [];
for i = 1:length(stats)
    if isempty(stats(i).name) == 1, continue; end
    n_sz = [n_sz;length(pt(i).sz)];
end

fprintf('The mean number of seizures was %1.1f (range %d-%d across patients).\n',...
    mean(n_sz),min(n_sz),max(n_sz));



end