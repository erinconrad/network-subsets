function network_stats(A)

%{

This function takes an adjacency matrix A, calculates the control
centrality of each node, and then randomly resamples the network, retaining
specified fractions of electrodes, and then gets statistics on how the
control centralities change.

%}

tic

%% Parameters
% How many random resamples to do of each fraction
n_perm = 1e2;

% What fraction of nodes to retain
e_f = [0.2 0.4 0.6 0.8 1];
n_f = length(e_f);

if isempty(A) == 1
    fprintf('A empty, using fake data.\n');
    [A,~] = makeFakeA(80,0.1);
    fprintf('%1.1f%% of possible links are connected.\n',...
        sum(sum(A==1))/size(A,1)^2*100);
end



%% Get true control centrality
c_c = control_centrality(A);
fprintf('There are %d synchronizing and %d desynchronizing nodes.\n',...
    sum(c_c<0),sum(c_c>0));

%% Resample network and get control centralities
% all_c_c is nch x n_f x n_perm size matrix
all_c_c = resampleNetwork(A,n_perm,e_f);

%% Initialize SMC and rho arrays
SMC = zeros(n_f,n_perm);
rho = zeros(n_f,n_perm);
true_cc_most_sync = zeros(n_f,n_perm);

%% Loop over each fraction and get various stats
for f = 1:n_f
    c_c_f = squeeze(all_c_c(:,f,:));
    
    for i_p = 1:n_perm
        c_c_f_p = squeeze(c_c_f(:,i_p));
        
        %% Do stats on this permutation of c_c and real c_c
        
        % if c_c_f_p all nans, skip it
        if sum(isnan(c_c_f_p)) == length(c_c_f_p)
            SMC(f,i_p) = nan;
            rho(f,i_p) = nan;
            continue
        end
        
        
        % Spearman rank coefficient
        rho_temp = corr(c_c(~isnan(c_c_f_p)),c_c_f_p(~isnan(c_c_f_p)),...
            'Type','Spearman');
        
        % Get binary versions for Simple matching coefficient
        c_c_bin = c_c(~isnan(c_c_f_p)) > 0;
        c_c_f_p_bin = c_c_f_p(~isnan(c_c_f_p)) > 0;
        
        % Get simple matching coefficient
        SMC_temp = simple_matching_coefficient(c_c_bin,c_c_f_p_bin);
        
        % Get the identity of the most synchronizing node (the one we would
        % tell the surgeons to resect)
        [~,ch_most_sync] = min(c_c_f_p);
        
        % Fill up SMC and rho arrays
        SMC(f,i_p) = SMC_temp;
        rho(f,i_p) = rho_temp;
        true_cc_most_sync(f,i_p) = c_c(ch_most_sync);
    end
    
end

%% Get mean and SD for SMC and rho for each fraction
SMC_mean = nanmean(SMC,2);
SMC_std = nanstd(SMC,0,2);

rho_mean = nanmean(rho,2);
rho_std = nanstd(rho,0,2);

%% How often would we resect the wrong piece of brain?
% For each f, get the % of times most synchronizing node is actually
% desynchronizing in the original network (how often would we tell the
% surgeon to resect a piece of brain that is thought to be desynchronizing
% in the original network)

resect_wrong = sum((true_cc_most_sync > 0),2)/n_perm;

%% Plot rho and SMC mean and std for each fraction
figure
set(gcf,'Position',[50 389 1400 409]);
subplot(1,3,1)
errorbar(e_f,rho_mean,rho_std,'linewidth',2);
xlabel('Fraction of original network included');
ylabel('Spearman rank coefficient');
title({'Spearman rank coefficient between original CC',...
    'and updated CC as a function of fraction of original network included'});

subplot(1,3,2)
errorbar(e_f,SMC_mean,SMC_std,'linewidth',2);
xlabel('Fraction of original network included');
ylabel('Simple matching coefficient');
title({'Simple matching coefficient between original CC',...
    'and updated CC as a function of fraction of original network included'});

subplot(1,3,3)
plot(e_f,resect_wrong*100,'linewidth',2);
xlabel('Fraction of original network included');
ylabel('% of permutations');
title({'Percent of time a desynchronizing node is',...
    'labeled as the most synchronizing'});

toc

end