function [adj,out_locs] = reconcileAdj(pt,whichPt)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
    pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;


 % Get name
name = pt(whichPt).name;

%% Load adjacency matrix
baseFolder = [mainFolder,'/data/adjacencyMatrices/',name,'/'];
listing = dir([baseFolder,'*.mat']);
load([baseFolder,listing.name]);
elecs = adj(7).data;

% Get the names of the unignored channels, should be in the right order
out_names = elecs.labels(elecs.ignore == 0)';

%% Get locs
locs = pt(whichPt).electrodeData.locs;
names = pt(whichPt).electrodeData.names;
out_locs = [];

for i = 1:size(locs,1)
    if ismember(names{i,1},out_names) == 1
        out_locs = [out_locs;locs(i,1:3)];
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