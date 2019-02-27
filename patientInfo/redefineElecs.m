function pt = redefineElecs

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Load pt struct
load([dataFolder,'structs/info.mat']);



for whichPt = 1:length(pt)
    
    name = pt(whichPt).name;
    
    % Load adjacency matrix
    baseFolder = [mainFolder,'/data/adjacencyMatrices/',name,'/'];
    listing = dir([baseFolder,'*.mat']);
    load([baseFolder,listing.name]);
    elecs = adj(7).data;
    
    % Setup new electrodeData
    pt(whichPt).new_elecs = [];
    base = pt(whichPt).new_elecs;
    base.locs = [];
    base.names = {};
    adj_names = {};
    count = 0;
    
    % Loop through electrodes in adjacency matrix and populate them into
    % new electrode data IF they are not ignored in adjacency matrix
    for i = 1:length(elecs.ignore)
        if elecs.ignore(i) == 0
            count = count + 1;
            base.electrodes(count).name = elecs.labels{i};
            adj_names = [adj_names;elecs.labels{i}];
            
            % Find the matching electrode in the original electrode struct
            % and get info about it
            found_it = 0;
            for j = 1:length(pt(whichPt).electrodeData.electrodes)
                if strcmp(pt(whichPt).electrodeData.electrodes(j).name,elecs.labels{i}) == 1
                    found_it = 1;
                    base.electrodes(count).loc = pt(whichPt).electrodeData.electrodes(j).xyz;
                    base.electrodes(count).type = pt(whichPt).electrodeData.electrodes(j).type;
                    base.locs = [base.locs;base.electrodes(count).loc];
                    base.names = [base.names;base.electrodes(count).name];
                    break
                end
            end
            
            % Throw an error if we couldn't find the electrode
            if found_it == 0
                error('Did not find %s\n', elecs.labels{i});
            end
            
        end
    end
    
    
    %% Confirm the order of electrodes is identical
    if isequal(adj_names,base.names) == 0
        error('Names don''t line up!\n');
    end
    
    %% Open up electrode data file and confirm locs match
    electrodeFile = pt(whichPt).electrode_labels;
    if isempty(electrodeFile) ==1
        error('No electrode data!\n');
    end
    
    fileID = fopen([electrodeFolder,electrodeFile]);

    out=textscan(fileID, '%s', 'whitespace',',');
    out = out{1};
    nChans = length(out)/8;
    
    % Loop through the electrodes in the electrode file
    for k = 1:nChans
        
        % Get the name and location of the electrode
        temp_name = out{(k-1)*8+5};
        temp_loc = [str2double(out{(k-1)*8+2}) str2double(out{(k-1)*8+3}) str2double(out{(k-1)*8+4})];
        found_it = 0;
        
        % Loop through my new electrode struct
        for j = 1:length(base.electrodes)
            
            if strcmp(temp_name,base.electrodes(j).name) == 1
                found_it = 1;
                
                % confirm that locs match
                if isequal(temp_loc,base.electrodes(j).loc) == 0
                    error('Locations don''t match!\n');
                end
                
                break
                
            end
                
        end
        
        if found_it == 0
            % if we cannot find the electrode, confirm that it is an
            % ignored electrode in the adjacency matrix, or not in the
            % adjacency matrix
            
            found_in_adj = 0;
            for j = 1:length(elecs.labels)
                if strcmp(elecs.labels{j},temp_name) == 1
                    found_in_adj = 1;
                    if elecs.ignore(j) == 0
                        error('Electrode %s not ignored!\n',temp_name);
                    end
                    break
                end
                
            end
            if found_in_adj == 0
                fprintf('Did not include %s as not in adj\n',temp_name);
            end
        end
        
    end
    close(fid)
    
    pt(whichPt).new_elecs = base;
end

end