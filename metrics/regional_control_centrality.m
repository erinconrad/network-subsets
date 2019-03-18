function [cc,elecs] = regional_control_centrality(A,n,locs,contig)

%{
This function calculates control centrality of a region rather than a
single node. 

Input parameters:
A - the adjacency matrix
n - the number of electrodes in the region
locs - the locations of all electrodes
contig - 1 if contiguous electrodes, 0 if random


Output parameters:
cc - an array with the control centrality of each region
elecs -  an array with the identities of electrodes in this region

Note: There are N_total_electrodes choose n_region_electrodes possible
choices of regions. If N_total_electrodes is 100 and n_region_electrodes is
20, then this means there are 5x10^20 ways to choose regions. This is too
many to enumerate all possible regions and find the lowest possible cc. 

I am going to reduce 5x10^20 to just N_total_electrodes by looping
through all electrodes, calculating its nearest n_region_electrodes-1
neighbors, and getting control centrality for that region.

Other future options would be to try to enumerate all possible contiguous
subsets (not sure how many of these there are, probably a lot), or to do a
Monte Carlo approach and randomly choose ~10^5 subsets of electrodes
(perhaps with spatial constraints) and pick whichever one is lowest from
that Monte Carlo approach and say that has the lowest control centrality.

%}

% Calculate the synchronizability
sync = synchronizability(A);

if contig == 0
    
    %{
    % I could do Monte Carlo approach here. Warning! If I do 1e4, I get a
    % different set of "min cc electrodes" each time. 
    nboot = 1e4;
    
    cc = nan(nboot,1);
    elecs = nan(nboot,n);
    
    for ib = 1:nboot
        if mod(ib,1e3) == 0
            fprintf('%d\n',ib);
        end
        nchs = size(locs,1);
        chs = 1:nchs;
        %% Build random set of "contiguous" electrodes
        % We will start with a random electrode and add its nearest
        % electrode. Then we will pick a random one of those 2 and add its
        % nearest. Then we will pick a random one of those 3 and add its
        % nearest, etc.
        
        % Get first electrode
        set = randi(nchs);
        
        while 1
            % pick a random electrode in the set
            curr_elec = set(randi(length(set)));
            
            % Get distances to other chs
            dist = vecnorm(locs-repmat(locs(curr_elec,:),nchs,1),2,2);
            [~,I] = sort(dist);
            chs_sorted = chs(I);

            % loop through these and add the first one that is not already
            % in the set to the set
            for i = 1:length(chs_sorted)
                if ismember(chs_sorted(i),set) == 0
                    set = [set;chs_sorted(i)];
                    break
                end
            end
           
            % break if reached size
            if length(set) == n
                break
            end
        end
        
        if 1 == 0
            % plot the set
            figure
            scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
            hold on
            scatter3(locs(set,1),locs(set,2),locs(set,3),100,'r','filled');
            pause
            close(gcf)
        end
            
        % Remove the channels
        A_temp = A;
        A_temp(set,:) = [];
        A_temp(:,set) = [];
        
        % Recalculate the synchronizability
        sync_temp = synchronizability(A_temp);
        
        % Calculate control centrality
        cc(ib) = (sync_temp - sync)/sync;
        
        elecs(ib,:) = set';
    end
    %}
    
elseif contig == 1
    % number of total channels
    nchs = size(locs,1);
    chs = 1:nchs;
    
    % initialize arrays
    elecs = nan(nchs,n);
    cc = nan(nchs,1);
    
    for i = 1:nchs
        
        % Get distances to other chs
        dist = vecnorm(locs-repmat(locs(i,:),nchs,1),2,2);
        
        % sort channels by distahces
        [~,I] = sort(dist);
        chs_sorted = chs(I);
        
        % Take the first n channels
        region = chs_sorted(1:n);
        elecs(i,:) = region;
        
        if 1 == 0
            % plot the set
            figure
            scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
            hold on
            scatter3(locs(region',1),locs(region',2),locs(region',3),100,'r','filled');
            pause
            close(gcf)
        end
        
        % Remove the channels
        A_temp = A;
        A_temp(region,:) = [];
        A_temp(:,region) = [];
        
        % Recalculate the synchronizability
        sync_temp = synchronizability(A_temp);
        
        % Calculate control centrality
        cc(i) = (sync_temp - sync)/sync;
        
        
    end
    
    % Plot the regional control centralities
    if 1 == 0
        figure
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
        hold on
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,cc,'filled');
        colorbar;
        pause
        close(gcf)
    end
end


end