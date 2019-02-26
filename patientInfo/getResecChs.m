function pt = getResecChs(pt)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
    pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;

resecFolder = [electrodeFolder,'resectedElecs/for_erin/'];

for whichPt = 1:length(pt)
    name = pt(whichPt).name;
    
    % Find appropriate file in resected folder
    listing = dir([resecFolder,name,'*']);
    
    if length(listing) ~= 1
        fprintf('Warning, no resected electrode file for %s, skipping...\n',...
            name);
        continue
    end
    
    % Open file
    T = readtable([resecFolder,listing.name],'Delimiter',',',...
        'ReadVariableNames',false);
    
    new_elecs = [];
    new_labels = {};
    
    for i = 1:size(T,1)
        
        % Get the label and add it to my cell array
        elec_label = T.Var2(i);
        new_labels = [new_labels;elec_label];
        
        found_it = 0;
        
        % Find the corresponding electrode number in my electrode field
        for j = 1:length(pt(whichPt).electrodeData.electrodes)
            
            
            if strcmp(elec_label,pt(whichPt).electrodeData.electrodes(j).name) == 1
                % Add the electrode number to my array
                new_elecs = [new_elecs;j];
                found_it = 1;
                
                % break out of inner loop
                break
                
            end

        end
        
        if found_it == 0
            error('Did not find electrode!\n');
        end
        
    end
    
    %% Add them to my structure
    pt(whichPt).resec.names = new_labels;
    pt(whichPt).resec.nums = new_elecs;
    
    %% Confirm that they line up
    for j = 1:length(pt(whichPt).resec.nums)
        which_num = pt(whichPt).resec.nums(j);
        if strcmp(pt(whichPt).electrodeData.electrodes(which_num).name,pt(whichPt).resec.names{j}) == 0
            error('Electrode names and nums do not line up for %s\n',...
                pt(whichPt).name);
        end
    end
    
end
    
end

