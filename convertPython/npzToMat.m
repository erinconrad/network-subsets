function npzToMat(pt,whichPts)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
    pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
baseFolder = [mainFolder,'/data/adjacencyMatrices/'];

if isempty(whichPts) == 1
    whichPts = 1:length(pt);
end

for whichPt = whichPts
    
    % Get name
    name = pt(whichPt).name;
    
    % Get folder
    outputFolder = [baseFolder,name,'/'];
    
    if exist(outputFolder,'dir') == 0
        continue
    end
    
    if exist([outputFolder,'adj.mat'],'file') == 1
        continue
    end
    
    % Find the file that ends in .npz
    listing = dir([outputFolder,'*.npz']);
    fname = listing.name;
    
    % Unzip it
    if exist([outputFolder,'labels.npy'],'file') == 0
        unzip([outputFolder,fname],outputFolder);
    end

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
    
    % Add python module
    %{
    if count(py.sys.path,'') == 0
        insert(py.sys.path,int32(0),'');
    end
    %}
    
    new_file = [outputFolder,'labels.npy'];
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
    
end

end
