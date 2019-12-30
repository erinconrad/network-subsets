function get_variance(whichPts,which_sz)

%{
The purpose of this function is to get the distribution across electrodes
of variance in resampling, and to see if the distribution of nodal
properties is normal
%}

%% Parameters
e_f = [0.8]; % resample removing 20% of electrodes
freq_cell = {'high_gamma'};
which_times = 0;
contigs = 0;

%% Initialize things
n_f = length(e_f);
if isempty(whichPts) == 1
    whichPts = 1:33;
end


%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

pt = load([dataFolder,'structs/info.mat']);
pt = pt.pt;

out_folder = [resultsFolder,'variances/'];
if exist(out_folder,'dir') == 0
    mkdir(out_folder)
end

%% Loop through patients, times, frequencies
for ff = 1:length(freq_cell)
    
    freq = freq_cell{ff};
    fprintf('Doing %s\n',freq);
    
    for which_sec = which_times
        fprintf('Doing %d second\n',which_sec);
        
        if which_sec < 0
            sec_text = sprintf('sec_neg%d',abs(which_sec));
        else
            sec_text = sprintf('sec_%d',abs(which_sec));
        end
        
        for whichPt = whichPts
            
            fprintf('Doing patient %d\n', whichPt);
            
            % Get locations
            locs = pt(whichPt).new_elecs.locs; % all electrode locations
            
            for contig = contigs
                fprintf('Doing contig %d\n',contig);
                % Skip if all electrode locations are -1 (means we don't have electrode
                % locations)
                if unique(pt(whichPt).new_elecs.locs) == -1
                    continue
                end
                
                if contig == 1
                    contig_text = 'contiguous';
                elseif contig == 0
                    contig_text = 'random';
                end
                
                if contig == 1
                    % Here, not taking random samples, but rather systematically going
                    % through each electrode and its N nearest neighbors
                    n_perm = length(pt(whichPt).new_elecs.electrodes);
                elseif contig == 0 || contig == 2 || contig == 4
                    % Take 1000 random permutations
                    n_perm = 1e3;
                end
                
                %% Get adjacency matrix and get right frequency and time
                [adj,~,sz_num] = reconcileAdj(pt,whichPt,which_sz);
                
                
                if isempty(adj) == 1
                    fprintf('Cannot do %s\n\n',name);
                    continue;
                end
                
                if strcmp(freq,'high_gamma') == 1
                    A_all = adj(4).data;
                    if contains(adj(4).name,'highgamma') == 0
                        error('This isn''t gamma!'\n');
                    end
                end
                
                % Start with the middle and add which second
                if ceil(size(A_all,1)/2)+which_sec <= 0, continue; end
                if ceil(size(A_all,1)/2)+which_sec > size(A_all,1), continue; end
                A = squeeze(A_all(ceil(size(A_all,1)/2)+which_sec,:,:));
                if sum(sum(isnan(A))) == sum(sum(ones(size(A))))
                    continue
                end
                
                %% Get true metrics
                c_c = control_centrality(A);
                ns = node_strength(A);
                bc = betweenness_centrality(A,1);
                ec = eigenvector_centrality_und(A);
                clust = clustering_coef_wu(A);
                sync = synchronizability(A);
                eff = efficiency_wei(A, 0);
                trans = transitivity_wu(A);
                
                out(whichPt).(freq).(contig_text).(sec_text).metric(1).name = 'control centrality';
                out(whichPt).(freq).(contig_text).(sec_text).metric(2).name = 'node strength';
                out(whichPt).(freq).(contig_text).(sec_text).metric(3).name = 'betweenness centrality';
                out(whichPt).(freq).(contig_text).(sec_text).metric(4).name = 'eigenvector centrality';
                out(whichPt).(freq).(contig_text).(sec_text).metric(5).name = 'clustering coefficient';
                out(whichPt).(freq).(contig_text).(sec_text).metric(6).name = 'synchronizability';
                out(whichPt).(freq).(contig_text).(sec_text).metric(7).name = 'global efficiency';
                out(whichPt).(freq).(contig_text).(sec_text).metric(8).name = 'transitivity';
                
                out(whichPt).(freq).(contig_text).(sec_text).metric(1).true = c_c;
                out(whichPt).(freq).(contig_text).(sec_text).metric(2).true = ns;
                out(whichPt).(freq).(contig_text).(sec_text).metric(3).true = bc;
                out(whichPt).(freq).(contig_text).(sec_text).metric(4).true = ec;
                out(whichPt).(freq).(contig_text).(sec_text).metric(5).true = clust;
                out(whichPt).(freq).(contig_text).(sec_text).metric(6).true = sync;
                out(whichPt).(freq).(contig_text).(sec_text).metric(7).true = eff;
                out(whichPt).(freq).(contig_text).(sec_text).metric(8).true = trans;
                
                %% Plot histogram
                % Appear to be fairly non normal on a patient specific
                % level
                if 0
                    for i = 1:length(out(whichPt).(freq).(contig_text).(sec_text).metric)
                        figure
                        histogram(out(whichPt).(freq).(contig_text).(sec_text).metric(i).true)
                        title(sprintf('%s',out(whichPt).(freq).(contig_text).(sec_text).metric(i).name))
                        pause
                        close(gcf)
                    end
                end
                
                
                %% Resample network and get new metrics
                % all_c_c is nch x n_f x n_perm size matrix
                [all_c_c,all_ns,all_bc,all_sync,all_eff,overlap_soz,dist_soz,...
                    overlap_resec,dist_resec,elecs_min,...
                    all_par,all_trans,avg_par_removed,avg_bc_removed,...
                    all_sync_norm,all_eff_norm,all_trans_norm,all_ec,...
                    all_clust,all_le,cc_reg,dist_nearest_resec,sz_soz_dist,...
                    all_cc_norm,all_ns_norm,all_bc_norm,all_ec_norm,all_clust_norm] = ...
                    resampleNetwork(A,n_perm,e_f,contig,pt,whichPt,adj,sz_num);
                
                %% Fill up 
                resamp(1).val = all_c_c;
                resamp(2).val = all_ns;
                resamp(3).val = all_bc;
                resamp(4).val = all_ec;
                resamp(5).val = all_clust;
                resamp(6).val = all_sync;
                resamp(7).val = all_eff;
                resamp(8).val = all_trans;
                
                %% Get variance across resamples (only makes sense for nodal metrics)
                for i = 1:5
                    var = squeeze(nanstd(resamp(i).val,0,3).^2);
                    
                    % Store it
                    out(whichPt).(freq).(contig_text).(sec_text).metric(i).var_resamp = var;
                end
                
                % Save
                save([out_folder,'variances.mat'],'out')
                
            end
            
        end
    end
end


end