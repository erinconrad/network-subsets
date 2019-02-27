function network_stats(whichPts)

%{

This function takes an adjacency matrix A, calculates the control
centrality of each node, and then randomly resamples the network, retaining
specified fractions of electrodes, and then gets statistics on how the
control centralities change.


Stuff to add:
- try a series of adjacency matrices, try to replicate the main finding of
the VCR paper with subsampling
- get other network metrics 
    - path length
    (https://www.sciencedirect.com/science/article/pii/S1388245715012584)
    - clustering coefficient
    (https://www.sciencedirect.com/science/article/pii/S1388245715012584)
    - node heterogeneity
    - epileptor model - The Virtual Brain

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

% 1 if doing SOZ analysis, 0 if doing main analysis
do_soz_analysis = 1;

doPlots = 1;

% add to existing stats array?
merge = 1;

% How many random resamples to do of each fraction
n_perm = 1e2;

% What fraction of nodes to retain
if do_soz_analysis == 1
    e_f = 0.8;
else
    e_f = [0.2 0.4 0.6 0.8 1];
end
n_f = length(e_f);

% Which freq?
freq = 'high_gamma';


%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

load([dataFolder,'structs/info.mat']);

if merge == 1
    if do_soz_analysis == 1
        if exist([resultsFolder,'basic_metrics/soz.mat'],'file') ~= 0
            load([resultsFolder,'basic_metrics/soz.mat']);
        else
            soz = struct;
        end
    else
        if exist([resultsFolder,'basic_metrics/stats.mat'],'file') ~= 0
            load([resultsFolder,'basic_metrics/stats.mat']);
        else
            stats = struct;
        end
    end
else
    if do_soz_analysis == 1
        soz = struct;
    else
        stats = struct;
    end
end


%% Loop through patients, times, and whether contig or random electrodes
for which_sec = [-5 0] % 0 means start time of the seizure, -5 is 5 seconds before
for contig = [1 0]

for whichPt = whichPts

    % Make result folder
    name = pt(whichPt).name;
    
    fprintf('Doing %s\n',name);
    
    outFolder = [resultsFolder,'basic_metrics/',name,'/'];
    if exist(outFolder,'dir') == 0
        mkdir(outFolder);
    end

    if contig == 1
        contig_text = 'contiguous';
    else
        contig_text = 'random';
    end
    
    if which_sec < 0
        sec_text = sprintf('sec_neg%d',abs(which_sec));
    else
        sec_text = sprintf('sec_%d',abs(which_sec));
    end
    
    % Continue if we've already done it
    if merge == 1
        if do_soz_analysis == 1
            stats = soz;
        end
        if length(stats) >= whichPt
            if isfield(stats(whichPt),'eff') == 1
                if isfield(stats(whichPt).eff,(contig_text)) == 1
                    if isfield(stats(whichPt).eff.(contig_text),(sec_text)) == 1
                        if isfield(stats(whichPt).eff.(contig_text).(sec_text),'true') == 1
                            if isempty(stats(whichPt).eff.(contig_text).(sec_text).true) == 0
                                fprintf('Did %s, skipping\n',name);
                              %  continue
                            end
                        end
                    end
                end
            end 
        end
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
    [adj,~] = reconcileAdj(pt,whichPt);

    if strcmp(freq,'high_gamma') == 1
        A_all = adj(4).data;
    end

    % Start with the middle and add which second
    A = squeeze(A_all(ceil(size(A_all,1)/2)+which_sec,:,:));

    %% Get true control centrality
    c_c = control_centrality(A);
    fprintf('There are %d synchronizing and %d desynchronizing nodes.\n',...
        sum(c_c<0),sum(c_c>0));
    
    % Get identity of node with lowest control centrality
    [~,min_cc_true] = min(c_c);
    
    % Get location of node with lowest control centrality
    locs = pt(whichPt).new_elecs.locs; % all electrode locations
    min_cc_true_loc = locs(min_cc_true,:);

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
    [all_c_c,all_ns,all_bc,all_sync,all_eff,overlap_soz,dist_soz,...
        overlap_resec,dist_resec] = ...
        resampleNetwork(A,n_perm,e_f,contig,pt,whichPt,adj);

    %% Initialize SMC and rho arrays for node-level metrics

    % Control centrality stuff
    SMC_cc = zeros(n_f,n_perm);
    rho_cc = zeros(n_f,n_perm);
    true_cc_most_sync = zeros(n_f,n_perm);
    dist_cc = zeros(n_f,n_perm);

    % Betweenness centrality stuff
    rho_bc = zeros(n_f,n_perm);

    % Node strength centrality stuff
    rho_ns = zeros(n_f,n_perm);

    %% Loop over each fraction and get various stats
    for f = 1:n_f
        c_c_f = squeeze(all_c_c(:,f,:));
        ns_f = squeeze(all_ns(:,f,:));
        bc_f = squeeze(all_bc(:,f,:));
        
        % Loop over each permutation
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
                
                % Get distances between lowest control centrality electrode
                % in resampled network and original network
                [~,min_cc_resample] = min(c_c_f_p);
                min_cc_resample_loc = locs(min_cc_resample,:);
                dist_cc(f,i_p) = sqrt(sum((min_cc_true_loc-min_cc_resample_loc).^2));

            end


            %% Do Spearman rank on node strength and betweenness centrality
            % For these, they are always non-negative and so SMC doesn't make
            % sense
            [rho_ns(f,i_p),~] = doStats(ns,ns_f_p);
            [rho_bc(f,i_p),~] = doStats(bc,bc_f_p);


        end

    end

    %% Get mean and SD for SMC and rho for each fraction
    % This uses Fisher's transformation for the Spearman rank coefficients
    % but nothing special for SMC
    
    SMC_mean_cc = nanmean(SMC_cc,2);
    SMC_std_cc = nanstd(SMC_cc,0,2);

    rho_mean_cc = average_rho(rho_cc,2);
    rho_std_cc = nanstd(rho_cc,0,2);

    rho_mean_ns = average_rho(rho_ns,2);
    rho_std_ns = nanstd(rho_ns,0,2);

    rho_mean_bc = average_rho(rho_bc,2);
    rho_std_bc = nanstd(rho_bc,0,2);

    %% How often would we resect the wrong piece of brain?
    % For each f, get the % of times most synchronizing node is actually
    % desynchronizing in the original network (how often would we tell the
    % surgeon to resect a piece of brain that is thought to be desynchronizing
    % in the original network)
    resect_wrong = sum((true_cc_most_sync > 0),2)/n_perm;
    
    
    %% Distance from truest min cc to min cc in resampled network
    cc_dist_mean = nanmean(dist_cc,2);
    cc_dist_std = nanstd(dist_cc,0,2);
    
    %% Fill up stats structure
    if do_soz_analysis == 0
    
        stats(whichPt).name = name;

        % control centrality
        stats(whichPt).cc.name = 'control centrality';
        stats(whichPt).cc.(contig_text).(sec_text).SMC.mean = SMC_mean_cc;
        stats(whichPt).cc.(contig_text).(sec_text).SMC.std = SMC_std_cc;
        stats(whichPt).cc.(contig_text).(sec_text).rho.mean = rho_mean_cc;
        stats(whichPt).cc.(contig_text).(sec_text).rho.std = rho_std_cc;
        stats(whichPt).cc.(contig_text).(sec_text).resect_wrong = resect_wrong;
        
        % Mean and STD of distance between true min cc and min cc in
        % resampled network
        stats(whichPt).cc.(contig_text).(sec_text).dist.mean = cc_dist_mean;
        stats(whichPt).cc.(contig_text).(sec_text).dist.std = cc_dist_std;
        

        % node strength
        stats(whichPt).ns.name = 'node strength';
        stats(whichPt).ns.(contig_text).(sec_text).rho.mean = rho_mean_ns;
        stats(whichPt).ns.(contig_text).(sec_text).rho.std = rho_std_ns;

        % betweenness centrality
        stats(whichPt).bc.name = 'betweenness centrality';
        stats(whichPt).bc.(contig_text).(sec_text).rho.mean = rho_mean_bc;
        stats(whichPt).bc.(contig_text).(sec_text).rho.std = rho_std_bc;

        % synchronizability
        stats(whichPt).sync.name = 'synchronizability';
        stats(whichPt).sync.(contig_text).(sec_text).mean = mean(all_sync,2);
        stats(whichPt).sync.(contig_text).(sec_text).std = std(all_sync,0,2);
        stats(whichPt).sync.(contig_text).(sec_text).true = sync;

        % efficiency
        stats(whichPt).eff.name = 'global efficiency';
        stats(whichPt).eff.(contig_text).(sec_text).mean = mean(all_eff,2);
        stats(whichPt).eff.(contig_text).(sec_text).std = std(all_eff,0,2);
        stats(whichPt).eff.(contig_text).(sec_text).true = eff;

        save([resultsFolder,'basic_metrics/stats.mat'],'stats');
    elseif do_soz_analysis == 1
        rel_sync = (all_sync-sync)/sync;
        rel_eff = (all_eff-eff)/eff;
        soz(whichPt).(contig_text).(sec_text).rho_cc = rho_cc;
        soz(whichPt).(contig_text).(sec_text).rho_bc = rho_bc;
        soz(whichPt).(contig_text).(sec_text).rho_ns = rho_ns;
        soz(whichPt).(contig_text).(sec_text).sync = rel_sync;
        soz(whichPt).(contig_text).(sec_text).eff = rel_eff;
        soz(whichPt).(contig_text).(sec_text).dist_soz = dist_soz;
        soz(whichPt).(contig_text).(sec_text).dist_resec = dist_resec;
        soz(whichPt).(contig_text).(sec_text).overlap_soz = overlap_soz;
        soz(whichPt).(contig_text).(sec_text).overlap_resec = overlap_resec;
        
        save([resultsFolder,'basic_metrics/soz.mat'],'soz');
    end

    %% Plots
  
    if doPlots == 1
        
        if do_soz_analysis == 1
        
        for i = 1:n_f
            %% Plot metrics as a function of relation to SOZ
            figure
            set(gcf,'Position',[92 126 1349 672]);
            % Distance from SOZ, global metrics
            subplot(2,2,1)
            scatter(dist_soz(i,:),rel_sync(i,:),100,'filled','b');
            hold on
            scatter(dist_soz(i,:),rel_eff(i,:),100,'filled','r');
         %   legend('Synchronizability','True synchronizability','Global efficiency','True global efficiency',...
         %       'location','northeastoutside');
            title('Global metrics as a function of distance from SOZ');
            xlabel('Distance from SOZ');
            ylabel('Metric');
            set(gca,'fontsize',20);

            % Distance from SOZ, nodal metrics
            subplot(2,2,2)
            scatter(dist_soz(i,:),rho_cc(i,:),100,'filled','b');
            hold on
            scatter(dist_soz(i,:),rho_bc(i,:),100,'filled','r');
            scatter(dist_soz(i,:),rho_ns(i,:),100,'filled','g');
            %legend('Control centrality','Betweenness centrality','Node strength','location','northeastoutside');
            title('Correlation of nodal metrics');
            xlabel('Distance from SOZ');
            ylabel('Correlation of metric');
            set(gca,'fontsize',20);

            % Overlap with SOZ, global metrics
            subplot(2,2,3)
            scatter(overlap_soz(i,:)*100,rel_sync(i,:),100,'filled','b');
            hold on
            scatter(overlap_soz(i,:)*100,rel_eff(i,:),100,'filled','r');
            %legend('Synchronizability','True synchronizability','Global efficiency','True global efficiency');
            title('Global metrics as a function of overlap with SOZ');
            xlabel('% of electrodes removed that were in SOZ');
            ylabel('Metric');
            set(gca,'fontsize',20);

            % Overlap with SOZ, nodal metrics
            subplot(2,2,4)
            scatter(overlap_soz(i,:)*100,rho_cc(i,:),100,'filled','b');
            hold on
            scatter(overlap_soz(i,:)*100,rho_bc(i,:),100,'filled','r');
            scatter(overlap_soz(i,:)*100,rho_ns(i,:),100,'filled','g');
            %legend('Control centrality','Betweenness centrality','Node strength');
            title('Correlation of nodal metrics');
            xlabel('% of electrodes removed that were in SOZ');
            ylabel('Correlation of metric');
            set(gca,'fontsize',20);
            print(gcf,[outFolder,'soz_',contig_text,sec_text],'-depsc');
            close(gcf)

            %% Plot metrics as a function of relation to resection
            figure
            set(gcf,'Position',[92 126 1349 672]);
            % Distance from resec, global metrics
            subplot(2,2,1)
            scatter(dist_resec(i,:),rel_sync(i,:),100,'filled','b');
            hold on
            scatter(dist_resec(i,:),rel_eff(i,:),100,'filled','r');
           % legend('Synchronizability','True synchronizability','Global efficiency','True global efficiency',...
           %     'location','northeastoutside');
            title('Global metrics');
            xlabel('Distance from resected electrodes');
            ylabel('Metric');
            set(gca,'fontsize',20);

            % Distance from resec, nodal metrics
            subplot(2,2,2)
            scatter(dist_resec(i,:),rho_cc(i,:),100,'filled','b');
            hold on
            scatter(dist_resec(i,:),rho_bc(i,:),100,'filled','r');
            scatter(dist_resec(i,:),rho_ns(i,:),100,'filled','g');
           % legend('Control centrality','Betweenness centrality','Node strength',...
           %     'location','northeastoutside');
            title('Correlation of nodal metrics');
            xlabel('Distance from resected electrodes');
            ylabel('Correlation of metric');
            set(gca,'fontsize',20);

            % Overlap with resec, global metrics
            subplot(2,2,3)
            scatter(overlap_resec(i,:)*100,rel_sync(i,:),100,'filled','b');
            hold on
            scatter(overlap_resec(i,:)*100,rel_eff(i,:),100,'filled','r');
           % legend('Synchronizability','True synchronizability','Global efficiency','True global efficiency');
            title('Global metrics');
            xlabel('% of electrodes removed that were in resected region');
            ylabel('Metric');
            set(gca,'fontsize',20);

            % Overlap with resec, nodal metrics
            subplot(2,2,4)
            scatter(overlap_resec(i,:)*100,rho_cc(i,:),100,'filled','b');
            hold on
            scatter(overlap_resec(i,:)*100,rho_bc(i,:),100,'filled','r');
            scatter(overlap_resec(i,:)*100,rho_ns(i,:),100,'filled','g');
           % legend('Control centrality','Betweenness centrality','Node strength');
            title('Correlation of nodal metrics');
            xlabel('% of electrodes removed that were in resected region');
            ylabel('Correlation');
            set(gca,'fontsize',20);
            print(gcf,[outFolder,'resec_',contig_text,sec_text],'-depsc');
            close(gcf)
        
        end
        else
        
        % Control centrality
        figure
        set(gcf,'Position',[50 389 1400 409]);
        subplot(1,3,1)
        errorbar(e_f,rho_mean_cc,rho_std_cc,'k','linewidth',2);
        xlabel('Fraction of original network included');
        ylabel('Spearman rank coefficient');
        title(sprintf(['Spearman rank coefficient between original CC\n'...
            'and updated CC as a function of fraction of original network included\n'...
            'taking %s electrodes, %s'],contig_text, sec_text));

        subplot(1,3,2)
        errorbar(e_f,SMC_mean_cc,SMC_std_cc,'k','linewidth',2);
        xlabel('Fraction of original network included');
        ylabel('Simple matching coefficient');
        title(sprintf(['Simple matching coefficient between original CC\n'...
            'and updated CC as a function of fraction of original network included\n'...
            'taking %s electrodes, %s'],contig_text,sec_text));

        subplot(1,3,3)
        plot(e_f,resect_wrong*100,'k','linewidth',2);
        xlabel('Fraction of original network included');
        ylabel('% of permutations');
        title(sprintf(['Percent of time a desynchronizing node is\n'...
            'labeled as the most synchronizing\n'...
            'taking %s electrodes, %s'],contig_text,sec_text));
        print(gcf,[outFolder,'cc_',contig_text,sec_text],'-depsc');
        close(gcf)


        % Node strength
        figure
        set(gcf,'Position',[50 389 500 409]);
        errorbar(e_f,rho_mean_ns,rho_std_ns,'k','linewidth',2);
        xlabel('Fraction of original network included');
        ylabel('Spearman rank coefficient');
        title(sprintf(['Spearman rank coefficient between original node strength\n'...
            'and updated node strength as a function of fraction of original network included\n'...
            'taking %s electrodes, %s'],contig_text,sec_text));
        print(gcf,[outFolder,'ns_',contig_text,sec_text],'-depsc');
        close(gcf)

        % Betweenness centrality
        figure
        set(gcf,'Position',[50 389 500 409]);
        errorbar(e_f,rho_mean_bc,rho_std_bc,'k','linewidth',2);
        xlabel('Fraction of original network included');
        ylabel('Spearman rank coefficient');
        title(sprintf(['Spearman rank coefficient between original BC\n'...
            'and updated BC as a function of fraction of original network included\n'...
            'taking %s electrodes, %s'],contig_text,sec_text));
        print(gcf,[outFolder,'bc_',contig_text,sec_text],'-depsc');
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
            'taking %s electrodes, %s'],contig_text,sec_text));
        print(gcf,[outFolder,'sync_',contig_text,sec_text],'-depsc');
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
            'taking %s electrodes, %s'],contig_text,sec_text));
        print(gcf,[outFolder,'eff_',contig_text,sec_text],'-depsc');
        close(gcf)
        end
    end
    
end

end
end


toc

end