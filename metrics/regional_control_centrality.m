function regional_control_centrality(A,n,locs,contig)

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

I am going to reduce 5x10^20 to just 100 (or N_total_electrodes) by looping
through all electrodes, calculating its nearest n_region_electrodes-1
neighbors, and getting control centrality for that region.

Other future options would be to try to enumerate all possible contiguous
subsets (not sure how many of these there are, probably a lot), or to do a
Monte Carlo approach and randomly choose ~10^5 subsets of electrodes
(perhaps with spatial constraints) and pick whichever one is lowest from
that Monte Carlo approach and say that has the lowest control centrality.

%}

if contig == 0
    error('I don''t know how to do this\n');
    % I could do Monte Carlo approach here.
    
elseif contig == 1
    
    %% Restrict to electrodes that appear in A
    % Some electrodes aren't in A because they don't have an EEG signal
    
    
end


end