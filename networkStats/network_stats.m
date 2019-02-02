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
        
        % Fill up SMC and rho arrays
        SMC(f,i_p) = SMC_temp;
        rho(f,i_p) = rho_temp;
        
    end
    
end

%% Get mean and SD for SMC and rho for each fraction
SMC_mean = nanmean(SMC,2);
SMC_std = nanstd(SMC,0,2);

rho_mean = nanmean(rho,2);
rho_std = nanstd(rho,0,2);

%% Plot rho and SMC mean and std for each fraction
figure
set(gcf,'Position',[440 389 1001 409]);
subplot(1,2,1)
errorbar(e_f,rho_mean,rho_std,'linewidth',2);
xlabel('Fraction of original network included');
ylabel('Spearman rank coefficient');
title({'Spearman rank coefficient between original CC',...
    'and updated CC as a function of fraction of original network included'});

subplot(1,2,2)
errorbar(e_f,SMC_mean,SMC_std,'linewidth',2);
xlabel('Fraction of original network included');
ylabel('Simple matching coefficient');
title({'Simple matching coefficient between original CC',...
    'and updated CC as a function of fraction of original network included'});

toc

end