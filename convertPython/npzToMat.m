function npzToMat(pt,whichPts)

mod = py.importlib.import_module('open_adj');
py.reload(mod);

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
    pwfile,dataFolder,bctFolder,mainFolder,adjFolder] = resectFileLocs;
baseFolder = [mainFolder,'/data/adjacencyMatrices/'];

if isempty(whichPts) == 1
    whichPts = 1:length(pt);
end

for whichPt = whichPts
    
    % Get name
    name = pt(whichPt).name;
    
    fprintf('Doing %s\n',name);
    
    % Get folder
    outputFolder = [baseFolder,name,'/'];
    
    % Current adjacency folder
    adj_pt_folder = [adjFolder,name,'/aim3/'];
    
    if exist(adj_pt_folder,'dir') == 0
        continue
    end
    
    if exist([outputFolder,'adj.mat'],'file') == 1
        continue
    end
    
    % Find the files that end in .npz
    listing = dir([adj_pt_folder,'/*multiband.npz']);
    
    % Find smallest number
    minNum = 1000;
    for n = 1:length(listing)
        fname = listing(n).name;
        [starti,endi] = regexp(fname,'Ictal.\d+.');
        which_multiband = fname(starti + 6:endi-1);
        which_mb_num = str2double(which_multiband);
        if which_mb_num < minNum
            whichFile = n;
            minNum = which_mb_num;
            which_mb_out = which_multiband;
        end
    end
    
    % Now find second smallest number
    smallest = minNum;
    minNum = 1000;
    which_mb_num = 1000;
    
    for n = 1:length(listing)
        fname = listing(n).name;
        [starti,endi] = regexp(fname,'Ictal.\d+.');
        which_multiband = fname(starti + 6:endi-1);
        which_mb_num = str2double(which_multiband);
        if which_mb_num < minNum && which_mb_num > smallest
            whichFile = n;
            minNum = which_mb_num;
            which_mb_out = which_multiband;
        end
    end
    
    if minNum == 1000
        fprintf('Warning, only one seizure for %s\n\n',name);
        continue;
    else
       
        fprintf('Doing seizure %d from %s\n',minNum,name);
        
    end
    
    
    fname = listing(whichFile).name;
    fprintf('Doing file %s\n\n',fname);
    
    
    % Unzip it
  %  if exist([outputFolder,'labels.npy'],'file') == 0
        unzip([adj_pt_folder,fname],outputFolder);
  %  end

    % Now find all the .npy files
    listing = dir([outputFolder,'*.npy']);
    
    count = 0;
    for j = 1:6
        fname = listing(j).name;
        filename = [outputFolder,fname];
        
        % define varname
        varname = fname(1:end-4);
        
        if strcmp(varname,'epoch_length') == 1
            continue
        end
        count = count+1;
        
        % Get header info
        [arrayShape, dataType, fortranOrder, littleEndian, totalHeaderLength, npyVersion] = readNPYheader(filename);
        
        % Get the data
        if fortranOrder
            data = memmapfile(filename, 'Format', {dataType, arrayShape, 'd'}, 'Offset', totalHeaderLength);

        else
            % Note! In this case, the dimensions of the array will be transposed,
            % e.g. an AxBxCxD array becomes DxCxBxA. 
            data = memmapfile(filename, 'Format', {dataType, arrayShape(end:-1:1), 'd'}, 'Offset', totalHeaderLength);

        end
        
       
        
        % add it to structure
        adj(count).name = varname;
        adj(count).data = data.Data.d;
        
    end
    
    %% Fix for the labels
    new_file = [outputFolder,'labels.npy'];
    label_data = py.open_adj.open_ad_f(new_file);
    vals = cell(label_data{1});
    vals = cell2mat(vals)';
    keys = cell(label_data{2});
    cellP = cell(1, numel(keys));
    for n= 1:numel(keys)
        strP = char(keys{n});
        cellP(n) = {strP};
    end
    keys = (cellP);
    
     % Re-sort by number
    [vals,I] = sort(vals);
    keys = keys(I);
    
    %% Figure out if we are ignoring it
    ignore = zeros(length(vals),1);
    for i = 1:length(ignore)
        if ismember(keys(i),pt(whichPt).ignore_electrodes) == 1
            ignore(i) = 1;
        end
    end
    
    % Add to structure
    adj(count+1).name = 'labels';
    adj(count+1).data.labels = keys;
    adj(count+1).data.nums = vals;
    adj(count+1).data.ignore = ignore;
    
    % Save the structure
    save([outputFolder,'adj',which_mb_out,'.mat'],'adj');
    

    % Add python module
    %{
    if count(py.sys.path,'') == 0
        insert(py.sys.path,int32(0),'');
    end
    %}
    
    %{
    
    label_data = char(py.str(py.open_adj.open_ad_f(new_file)));
    
    % Do regular expressions
    [start_num,end_num] = regexp(label_data,'(\d+,');
    [start_name,end_name] = regexp(label_data,'''\w*'':');
    ch_num = zeros(length(start_num),1);
    ch_name = cell(length(start_name),1);
    for k = 1:length(start_num)
        temp_num = str2double(label_data(start_num(k)+1:end_num(k)-1));
        ch_num(k) = temp_num;
        ch_name{k} = label_data(start_name(k)+1:end_name(k)-2);
    end
    
    % Add to structure
    adj(count+1).name = 'labels';
    adj(count+1).data.labels = ch_name;
    adj(count+1).data.nums = ch_num;
    
    % Save the structure
    save([outputFolder,'adj.mat'],'adj');
    %}
    
end

end
