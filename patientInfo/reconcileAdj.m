function [adj,out_locs] = reconcileAdj(pt,whichPt,which_sz)

%{
This function gets an adjacency matrix and gets the correct locations
corresponding with the row of each matrix
%}

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
    pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;


 % Get name
name = pt(whichPt).name;

%% Load adjacency matrix
baseFolder = [mainFolder,'/data/adjacencyMatrices/',name,'/'];
listing = dir([baseFolder,'*.mat']);


load([baseFolder,listing(which_sz).name]);
elecs = adj(7).data;

% Get the names of the unignored channels, should be in the same order as
% the rows in the adjacency matrix
out_names = elecs.labels(elecs.ignore == 0)';

%% Get locs
locs = pt(whichPt).new_elecs.locs;
names = pt(whichPt).new_elecs.names;
out_locs = [];

% Get the locations of the electrodes matching the names of the unignored
% electrodes in the adjacency matrix
for i = 1:size(locs,1)
    if ismember(names{i,1},out_names) == 1
        out_locs = [out_locs;locs(i,1:3)];
        
    else
        
        error('Warning! Could not find electrode %s in adjacency matrix structure for %s\n',...
            names{i,1},name);
        
    end
    
end

if size(out_names,1) ~= size(adj(1).data,2)
    error('What\n');
end

%% Check
%{
table(char(out_names))
table(char(elecs.labels'),elecs.ignore)
%}


end