function [adj,out_locs,sz_num] = reconcileAdj(pt,whichPt,which_sz,which_window)

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

if length(listing) < which_sz
    fprintf('Cannot do %s\n\n',name);
    adj = []; out_locs = []; sz_num = [];
    return
end

% get appropriate window
new_listing_count = 1;
for i = 1:length(listing)
    if contains(listing(i).name,'window') == 1
        if which_window ~=1
            new_listing(new_listing_count).name = listing(i).name;
            new_listing_count = new_listing_count + 1;
        end
    else
        if which_window == 1
            new_listing(new_listing_count).name = listing(i).name;
            new_listing_count = new_listing_count + 1;
        end
    end
end

% Confirm I don't have any "window" files when I shouldn't, and that all
% files have "window" if they should
for i = 1:length(new_listing)
    if which_window == 1
        if contains(new_listing(i).name,'window') == 1, error('what\n'); end
    else
        if contains(new_listing(i).name,'window') == 0, error('what\n'); end
    end
end

% note that "which_sz" is 1, 2, or 3 and whereas
% the sz number attached to the adjacency matrix is the actual seizure
% number. I need to make the conversion here
if which_window == 1
    load([baseFolder,new_listing(which_sz).name]);
elseif which_window == 2
    if contains(new_listing(1).name,'window2') == 0, error('what\n'); end
    load([baseFolder,new_listing(1).name]);
elseif which_window == 500
    if contains(new_listing(1).name,'window500') == 0, error('what\n'); end
    load([baseFolder,new_listing(500).name]);
end


elecs = adj(7).data;
sz_num = which_sz;
%{
s = regexp(listing(which_sz).name,'\d');

if whichPt == 20
    sz_num = which_sz;
    fprintf('Since %s, doing sz %d.\n',name,sz_num);
else
    sz_num = str2num(listing(which_sz).name(s));
    fprintf('This is seizure %d.\n',sz_num);
end
%}



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

if length(out_names) ~= size(adj(1).data,2)
    error('What\n');
end

%% Check
%{
table(char(out_names))
table(char(elecs.labels'),elecs.ignore)
%}


end