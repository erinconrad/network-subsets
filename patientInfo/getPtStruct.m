function getPtStruct

%% Parameters
outputFile = 'info.mat';

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,pwfile,dataFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);
ptInfo = loadjson(jsonfile);
ptnames = fieldnames(ptInfo.PATIENTS);

for i = 1:length(ptnames)

    info = ptInfo.PATIENTS.(ptnames{i});
    
    % Get basic info
    pt(i).name = ptnames{i};
    [pt(i).ieeg_name,electrodeFile,~,~,~] = ieegAndElectodeNames(pt(i).name);
    pt(i).chLocationFile = [pt(i).name,'_chLocations.mat'];
    pt(i).sz_onset = info.SeizureOnset;
    pt(i).clinical.outcome = info.Outcome;
    pt(i).clinical.lesionstatus = info.Lesion_0x20_Status;
    pt(i).clinical.pathology = info.Pathology;
    pt(i).clinical.resectiontype = info.Resection_0x20_Type;
    pt(i).clinical.ageOnset = info.AgeOnset;
    pt(i).clinical.ageSurgery = info.AgeSurgery;
    pt(i).clinical.sex = info.Sex;
    pt(i).clinical.seizureOnset = info.SeizureOnset;
    if isfield(info,'GridCoverageOfResection') == 1
        pt(i).clinical.gridcovresection = info.GridCoverageOfResection;
    end
    
    % Get electrode label file names
    pt(i).electrode_labels = electrodeFile;
    
    % Get seizures
    szs = fieldnames(info.Events.Ictal);
    
    whichSz = 0;
    
    for j = 1:length(szs)
        

       sz = info.Events.Ictal.(szs{j});
       
       if isfield(sz,'SeizureEEC') == 0
           continue
       end
       
       skipSz = 0;
       
   
       if sz.SeizureEEC == -1 || sz.SeizureEnd == -1
           skipSz = 1;
       end
       
       % If the seizure is the same time as a prior, skip this one
       if j > 1
           for k = 1:whichSz-1
               if abs(pt(i).sz(k).onset - sz.SeizureEEC) < 1
                  skipSz = 1;
                   
               end
               
           end
           
       end
       
      
       
       
       if skipSz == 1
           continue
       end
       
       whichSz = whichSz + 1;
       
       % Get seizure onset and offset
       pt(i).sz(whichSz).onset = sz.SeizureEEC;
       pt(i).sz(whichSz).offset = sz.SeizureEnd;
       pt(i).sz(whichSz).electrodes = sz.SEIZURE_ONSET_ELECTRODES;
       
       
       % If the seizure onset is before the prior seizure onset, switch
       % positions
       if j > 1
           if pt(i).sz(whichSz).onset < pt(i).sz(whichSz-1).onset
               firstSz = [pt(i).sz(whichSz).onset pt(i).sz(whichSz).offset];
               secondSz = [pt(i).sz(whichSz-1).onset pt(i).sz(whichSz-1).offset];
               pt(i).sz(whichSz-1).onset = firstSz(1);
               pt(i).sz(whichSz-1).offset = firstSz(2);
               
               pt(i).sz(whichSz).onset = secondSz(1);
               pt(i).sz(whichSz).offset = secondSz(2);
               
           end
           
       end
       
    end
    
    
end

%% Loop through patients, run times
for i = 1:length(pt)
    fprintf('Doing patient %s\n',pt(i).name);
    dataName =  pt(i).ieeg_name;
    if isempty(dataName) == 1
        continue
    end
    
    electrodeFile = pt(i).electrode_labels;
    if isempty(electrodeFile) ==1
        continue
    end
    
    %% Load EEG data info
    % calling this with 0 and 0 means I will just get basic info like sampling
    % rate and channel labels
    data = getiEEGData(dataName,0,0,pwfile);  


    %% Ignore certain electrodes (EKG, etc.)

    chLabels = data.chLabels;
    chLabelsParsed = chLabels;
    nchan = length(chLabelsParsed);
    chIgnore = zeros(nchan,1);

    % Get the channels labels from ieeg and then parse them to be
    % easier to read and compare to my names
    for ch = 1:length(chLabels)
       chLabelsParsed{ch} = chParser(chLabels{ch}); 
    end
    
    % Ignore elecs if they are not in the elec file
    ignoreElectrodes = findChsToIgnore(pt,i,chLabelsParsed);
    
    foundIgnoredChs = zeros(length(ignoreElectrodes),1);

    % Find channels that are equal to my ignored channels
    for x = 1:length(chLabelsParsed)
        chName = chLabelsParsed{x};
        for y = 1:length(ignoreElectrodes)
            if strcmp(chName,ignoreElectrodes{y}) == 1
                chIgnore(x) = 1;
                foundIgnoredChs(y) = 1;
            end
        end
    end

    % Give a warning if I have remaining ignoreElectrodes that I could
    % not identify amongst the channels
    if sum(foundIgnoredChs) ~= length(ignoreElectrodes)
        unfoundIdx = find(foundIgnoredChs==0);
        unfound = ignoreElectrodes{unfoundIdx};
        fprintf('Warning, could not find ignored channel %s\n',unfound);
    end

    channels = find(chIgnore == 0);
    unignoredChLabels = chLabelsParsed(channels);
    nchan = length(channels);
    electrodeData = chanLocUseGdf(unignoredChLabels,[electrodeFolder,electrodeFile]);
    electrodeData.allLabels = data.chLabels;

    % Add fs and electode data to the patient structure
    pt(i).fs = data.fs;
    pt(i).electrodeData = electrodeData;
    pt(i).electrodeData.locs = pt(i).electrodeData.locs(:,2:4);
    pt(i).channels = channels;
    
    % Get sz onset zone electrodes
    pt(i).newSOZChs = [];
    for j = 1:length(pt(i).sz)
        chnames = pt(i).sz(j).electrodes;
        chnums = zeros(length(chnames),1);
        for ich = 1:length(chnames)
            
                [Lia,chIds] = ismember(chnames{ich},pt(i).electrodeData.unignoredChs);
                if Lia == 0
                    fprintf('Warning, could not find channel %s in the list of unignored channels for patient %s\n',...
                        chnames{ich},pt(i).name);
                    error('');

                end
                chnums(ich) = chIds;
           



        end
        pt(i).sz(j).chs = chnums;
        pt(i).newSOZChs = [pt(i).newSOZChs;chnums];
    end
     pt(i).newSOZChs = unique(pt(i).newSOZChs);
    
    
end

save([dataFolder,'structs/',outputFile],'pt');


end