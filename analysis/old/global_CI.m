function global_CI(stats)

%% Parameters
contig_text = 'random';
sec_text = 'sec_0';
freq = 'high_gamma';
all_global = {'sync','eff','trans'};
all_global_names = {'Synchronizability','Global efficiency','Transitivity'};

%% Locations

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'cc_comparison/'];

figure
[ha,pos] = tight_subplot(2,1,[0.02 0.02],[0.02 0.02],[0.06 0.02]);

rho_perm = zeros(3,100);

for j = 1:length(all_global)
    
    
    true = zeros(29,1);
    all_perm = zeros(29,100);
    which_global = all_global{j};
    for i = 1:length(stats)    
        if isempty(stats(i).name) ==1
            true(i) = nan;
            all_perm(i,:) = nan(1,100);
            continue; 
        end

        base = stats(i).(freq).(contig_text).(sec_text);

        true(i) = base.(which_global).true;
        all_perm(i,:) = base.(which_global).all;

    end
    for i = 1:length(rho_perm)
        rho_perm(j,i) = corr(true(~isnan(true)),all_perm(~isnan(true),i),'Type','Spearman');
    end
    inan = find(isnan(true));
    true(inan) = [];
    all_perm(inan,:) = [];
    
    if j == 1
        axes(ha(1))
        violin(all_perm',true,'facecolor',[0 0.4470 0.7410])
        xticklabels([])

        ylabel(sprintf('%s',all_global_names{j}));


        xlabel('Which patient');
    end
    
end
%{
axes(ha(2))
for i = 1:size(rho_perm,1)
    scatter(i*ones(size(rho_perm,
end
%}



