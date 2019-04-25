function check_alignment(pt,whichPts)

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

for whichPt = whichPts
    name = pt(whichPt).name;
    baseFolder = [mainFolder,'/data/adjacencyMatrices/',name,'/'];
    listing = dir([baseFolder,'*.mat']);
    
    for i = 1:length(listing)
        s = str2num(regexp(listing(i).name,'\d'));
        if s ~= pt(whichPt).sz(i).adj_num
            error('what\n');
        end
        
        fprintf('Seizure %d aligned for %s\n',s,name);
    end
    
end


end