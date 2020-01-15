function unzip_npz_files(pt,whichPts,which_sz,window)


mod = py.importlib.import_module('open_adj');
py.reload(mod);

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
    pwfile,dataFolder,bctFolder,mainFolder,adjFolder,newAdjFolder] = resectFileLocs;
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
    if window == 1
        adj_pt_folder = [adjFolder,name,'/aim3/'];
    else
        adj_pt_folder = [newAdjFolder,name,'/'];
    end
    
    if exist(adj_pt_folder,'dir') == 0
        continue
    end
    
    if exist([outputFolder,'adj.mat'],'file') == 1
        continue
    end
    
    % Find the files that end in .npz
    if window == 1
        listing = dir([adj_pt_folder,'/*multiband.npz']);
    elseif window == 2
        listing = dir([adj_pt_folder,'/*multiband.2.npz']);
    elseif window == 500
        listing = dir([adj_pt_folder,'/*multiband.500.npz']);
    end
    
    if which_sz == 1
        % code to find the first seizure
        
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
        
        fprintf('Doing seizure %d from %s\n',minNum,name);
    
    elseif which_sz == 2
        
        % I believe this code will find the 2nd seizure if it exists
    
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
        end
        
        fprintf('Doing seizure %d from %s\n',minNum,name);
        
    elseif which_sz == 3 % this is code to do the last seizure
        
        if whichPt == 19
            
            % Fix - something wrong with last seizure; wrong # of
            % electrodes. I think it is because seizure #11 is the first sz
            % with the reimplanation. Seizure #6 is the last sz of the
            % original implant.
            name = 'HUP111A.Ictal.6.multiband.npz';
            which_mb_out = '6';
        else
        
        % Find the last seizure (that isn't 1 or 2)
        interictal_num = 1000;
        max_num = 1;
        which_mb_out = nan;
        which_mb_num = nan;
        
        for n = 1:length(listing)
            fname = listing(n).name;
            [starti,endi] = regexp(fname,'Ictal.\d+.');
            which_multiband = fname(starti + 6:endi-1);
            which_mb_num = str2double(which_multiband);
            
            if which_mb_num > max_num && which_mb_num < interictal_num
                whichFile = n;
                max_num = which_mb_num;
                which_mb_out = which_multiband;
            end
        end
        
        if isnan(which_mb_num) == 1 || exist([outputFolder,'adj',which_mb_out,'.mat'],'file') == 1 || which_mb_num == 1000
            fprintf('Warning, only one or two seizures for %s\n\n',name);
            continue
        end
        
        fprintf('Doing seizure %d from %s\n',max_num,name);
        
        end
        
    end
       
    
        
    
    if whichPt == 19 && last_sz == 1
        fname = name;
    else
        fname = listing(whichFile).name;
        fprintf('Doing file %s\n\n',fname);
    end

    
    
    % Unzip it
  %  if exist([outputFolder,'labels.npy'],'file') == 0
        unzip([adj_pt_folder,fname],outputFolder);
  %  end
  
    fprintf('Unzipped file %s to %s\n\n',fname,outputFolder);

    
    %% Loop through each npy file and make it a .mat file   
    py_file = 'save_as_mat.py';
    commandStr = ['python ',sprintf('%s "%s"',py_file,outputFolder)];
    system(commandStr);
    
    %% Load the mat file
    all = load([outputFolder,'all.mat']);
    
    %% Prepare adjacency matrices
    adj(1).name = 'all_adj_alphatheta';
    adj(1).data = permute(all.alphatheta,[3,1,2]);
    adj(2).name = 'all_adj_beta';
    adj(2).data = permute(all.beta,[3,1,2]);
    adj(3).name = 'all_adj_broadband_CC';
    adj(3).data = permute(all.broadband,[3,1,2]);
    adj(4).name = 'all_adj_highgamma';
    adj(4).data = permute(all.highgamma,[3,1,2]);
    adj(5).name = 'all_adj_lowgamma';
    adj(5).data = permute(all.lowgamma,[3,1,2]);
    adj(6).name = 'all_adj_veryhigh';
    adj(6).data = permute(all.veryhigh,[3,1,2]);
    
     %% Fix for the labels
    vals = all.label_vals;
    keys = all.label_keys;
    
    keys = cellstr(keys);
    
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
    adj(7).name = 'labels';
    adj(7).data.labels = keys;
    adj(7).data.nums = vals;
    adj(7).data.ignore = ignore;
    
    % Save the structure
    save([outputFolder,'adj',which_mb_out,'_window',window,'.mat'],'adj');
    
    % Delete all.mat (it causes all sorts of issues if I keep it there)
    delete([outputFolder,'all.mat'])
    
end

end
