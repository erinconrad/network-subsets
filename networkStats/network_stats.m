function network_stats(whichPt)

%{

This function takes an adjacency matrix A, calculates the control
centrality of each node, and then randomly resamples the network, retaining
specified fractions of electrodes, and then gets statistics on how the
control centralities change.


Stuff to add:
- try a series of adjacency matrices, try to replicate the main finding of
the VCR paper with subsampling
- get other network metrics 
    - node strength - I believe if you have adjacency matrix A, the node
    strength of i is sum(A(i,:)) --> could do SRC
    - betweenness centrality  --> could do SRC
    - path length
    (https://www.sciencedirect.com/science/article/pii/S1388245715012584)
    - clustering coefficient
    (https://www.sciencedirect.com/science/article/pii/S1388245715012584)
    - node heterogeneity
    - epileptor model - The Virtual Brain
    - overall synchronizability
    - global efficiency

- take number N of resected electrodes and randomly move them around so
still N contiguous electrodes and recalculate the control centrality of
the "resected region"

- which adjacency matrices to use?
    - 5 seconds before
    - right at seizure onset
    - 10 seconds after
- aim3/results/patientID/aim3/multiband - minutes long, nchxnchxTxfband
- multiple freq bands
    - start with high gamma

- brain connectivity toolbox

%}

tic

%% Parameters

% How many random resamples to do of each fraction
n_perm = 1e2;

% Remove a contiguous chunk of electrodes?
contig = 1;
if contig == 1
    contig_text = 'contiguous';
else
    contig_text = 'random';
end

% What fraction of nodes to retain
e_f = [0.2 0.4 0.6 0.8 1];
n_f = length(e_f);

% Which freq?
freq = 'high_gamma';

% Which second
which_sec = 0; % 0 means start time of the seizure, -10 is 10 seconds before

%% Load stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);

% Make result folder
name = pt(whichPt).name;
outFolder = [resultsFolder,'basic_metrics/',name,'/'];
if exist(outFolder,'dir') == 0
    mkdir(outFolder);
end
    

%{
if isempty(A) == 1
    fprintf('A empty, using fake data.\n');
    [A,~] = makeFakeA(size(locs,1),0.1);
    fprintf('%1.1f%% of possible links are connected.\n',...
        sum(sum(A==1))/size(A,1)^2*100);
end
%}

% Get locs and adjacency
[adj,locs] = reconcileAdj(pt,whichPt);

if strcmp(freq,'high_gamma') == 1
    A_all = adj(4).data;
end


A = squeeze(A_all(size(A_all,1)/2+which_sec,:,:));

%% Get true control centrality
c_c = control_centrality(A);
fprintf('There are %d synchronizing and %d desynchronizing nodes.\n',...
    sum(c_c<0),sum(c_c>0));

%% Get true synchronizability
sync = synchronizability(A);

%% Get true betweenness centrality
bc = betweenness_centrality(A,1);

%% Get true node strength
ns = node_strength(A);

%% Get true global efficiency
eff = efficiency_wei(A, 0);

%% Resample network and get new metrics
% all_c_c is nch x n_f x n_perm size matrix
[all_c_c,all_ns,all_bc,all_sync,all_eff] = ...
    resampleNetwork(A,n_perm,e_f,contig,locs);

%% Initialize SMC and rho arrays

% Control centrality stuff
SMC_cc = zeros(n_f,n_perm);
rho_cc = zeros(n_f,n_perm);
true_cc_most_sync = zeros(n_f,n_perm);

% Betweenness centrality stuff
rho_bc = zeros(n_f,n_perm);

% Node strength centrality stuff
rho_ns = zeros(n_f,n_perm);

%% Loop over each fraction and get various stats
for f = 1:n_f
    c_c_f = squeeze(all_c_c(:,f,:));
    ns_f = squeeze(all_ns(:,f,:));
    bc_f = squeeze(all_bc(:,f,:));
    
    for i_p = 1:n_perm
        c_c_f_p = squeeze(c_c_f(:,i_p));
        ns_f_p = squeeze(ns_f(:,i_p));
        bc_f_p = squeeze(bc_f(:,i_p));
        
        %% Do stats on control centrality
        
        % Get rho and SMC
        [rho_cc(f,i_p),SMC_cc(f,i_p)] = doStats(c_c,c_c_f_p);

        if sum(isnan(c_c_f_p)) == length(c_c_f_p)
            true_cc_most_sync(f,i_p) = nan;
        else
            
            % Get the identity of the most synchronizing node (the one we would
            % tell the surgeons to resect)
            [~,ch_most_sync] = min(c_c_f_p);

            % Fill up SMC and rho arrays
            true_cc_most_sync(f,i_p) = c_c(ch_most_sync);
            
        end

        
        %% Do stats on node strength and betweenness centrality
        % For these, they are always non-negative and so SMC doesn't make
        % sense
        [rho_ns(f,i_p),~] = doStats(ns,ns_f_p);
        [rho_bc(f,i_p),~] = doStats(bc,bc_f_p);
            
  
    end
    
end

%% Get mean and SD for SMC and rho for each fraction
SMC_mean_cc = nanmean(SMC_cc,2);
SMC_std_cc = nanstd(SMC_cc,0,2);

rho_mean_cc = nanmean(rho_cc,2);
rho_std_cc = nanstd(rho_cc,0,2);

rho_mean_ns = nanmean(rho_ns,2);
rho_std_ns = nanstd(rho_ns,0,2);

rho_mean_bc = nanmean(rho_bc,2);
rho_std_bc = nanstd(rho_bc,0,2);

%% How often would we resect the wrong piece of brain?
% For each f, get the % of times most synchronizing node is actually
% desynchronizing in the original network (how often would we tell the
% surgeon to resect a piece of brain that is thought to be desynchronizing
% in the original network)

resect_wrong = sum((true_cc_most_sync > 0),2)/n_perm;

%% Plots

% Control centrality
figure
set(gcf,'Position',[50 389 1400 409]);
subplot(1,3,1)
errorbar(e_f,rho_mean_cc,rho_std_cc,'k','linewidth',2);
xlabel('Fraction of original network included');
ylabel('Spearman rank coefficient');
title(sprintf(['Spearman rank coefficient between original CC\n'...
    'and updated CC as a function of fraction of original network included\n'...
    'taking %s electrodes'],contig_text));

subplot(1,3,2)
errorbar(e_f,SMC_mean_cc,SMC_std_cc,'k','linewidth',2);
xlabel('Fraction of original network included');
ylabel('Simple matching coefficient');
title(sprintf(['Simple matching coefficient between original CC\n'...
    'and updated CC as a function of fraction of original network included\n'...
    'taking %s electrodes'],contig_text));

subplot(1,3,3)
plot(e_f,resect_wrong*100,'k','linewidth',2);
xlabel('Fraction of original network included');
ylabel('% of permutations');
title(sprintf(['Percent of time a desynchronizing node is\n'...
    'labeled as the most synchronizing\n'...
    'taking %s electrodes'],contig_text));
print(gcf,[outFolder,'cc_',contig_text],'-depsc');
close(gcf)


% Node strength
figure
set(gcf,'Position',[50 389 500 409]);
errorbar(e_f,rho_mean_ns,rho_std_ns,'k','linewidth',2);
xlabel('Fraction of original network included');
ylabel('Spearman rank coefficient');
title(sprintf(['Spearman rank coefficient between original node strength\n'...
    'and updated node strength as a function of fraction of original network included\n'...
    'taking %s electrodes'],contig_text));
print(gcf,[outFolder,'ns_',contig_text],'-depsc');
close(gcf)

% Betweenness centrality
figure
set(gcf,'Position',[50 389 500 409]);
errorbar(e_f,rho_mean_bc,rho_std_bc,'k','linewidth',2);
xlabel('Fraction of original network included');
ylabel('Spearman rank coefficient');
title(sprintf(['Spearman rank coefficient between original BC\n'...
    'and updated BC as a function of fraction of original network included\n'...
    'taking %s electrodes'],contig_text));
print(gcf,[outFolder,'bc_',contig_text],'-depsc');
close(gcf)

% Synchronizability
figure
set(gcf,'Position',[50 389 500 409]);
errorbar(e_f,mean(all_sync,2),std(all_sync,0,2),'k','linewidth',2);
hold on
plot(get(gca,'xlim'),[sync sync],'k--','linewidth',2);
xlabel('Fraction of original network included');
ylabel('Synchronizability');
title(sprintf(['Synchronizability'...
    ' as a function of fraction of original network included\n'...
    'taking %s electrodes'],contig_text));
print(gcf,[outFolder,'sync_',contig_text],'-depsc');
close(gcf)

% Efficiency
figure
set(gcf,'Position',[50 389 500 409]);
errorbar(e_f,mean(all_eff,2),std(all_eff,0,2),'k','linewidth',2);
hold on
plot(get(gca,'xlim'),[eff eff],'k--','linewidth',2);
xlabel('Fraction of original network included');
ylabel('Global efficiency');
title(sprintf(['Global efficiency'...
    ' as a function of fraction of original network included\n'...
    'taking %s electrodes'],contig_text));
print(gcf,[outFolder,'eff_',contig_text],'-depsc');
close(gcf)


toc

end