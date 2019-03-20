function cc_npz_to_mat(pt,whichPts)

mod = py.importlib.import_module('open_cc');
mod_shape = py.importlib.import_module('get_shape');
py.reload(mod);
py.reload(mod_shape);

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
    pwfile,dataFolder,bctFolder,mainFolder,adjFolder] = resectFileLocs;
baseFolder = [mainFolder,'/data/control_centralities/'];

if isempty(whichPts) == 1
    whichPts = 1:length(pt);
end

for whichPt = whichPts
    
    % Get name
    name = pt(whichPt).name;
    
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
    
    listing = dir([adj_pt_folder,'/*1.noderes.npz']);
    fname = listing(1).name;
    
    % Unzip it
    unzip([adj_pt_folder,fname],outputFolder);
    
    % Now find the high gamma files
    listing = dir([outputFolder,'*control_centrality_highgamma.npy']);
    
    %% Load in all of the things
    fname = listing(1).name;
    fname = [outputFolder,fname];
    shape = cell(py.get_shape.get_shape_f(fname));
    
    nchs = shape{1};
    ntimes = shape{2};
    
    cc = nan(ntimes,nchs);
    
    for i = 1:ntimes
        data = cell(py.open_cc.open_cc_f(fname,1));
        data = cell2mat(data);
        
        cc(i,:) = data;
    end
    
    %% Save the array
    save([outputFolder,'cc_true.mat'],'cc');
    
end

end