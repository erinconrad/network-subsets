function pt = confirm_soz(pt)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;

ptInfo = loadjson(jsonfile);

ptnames = fieldnames(ptInfo.PATIENTS);

for whichPt = 1:length(pt)
    if isempty(pt(whichPt).sz) == 1, continue; end
    
    % Find the patient in the json file that matches
    for i = 1:length(ptnames)
        info = ptInfo.PATIENTS.(ptnames{i});
        if strcmp(pt(whichPt).name,ptnames{i}) == 1
            
            % Get seizures
            szs = fieldnames(info.Events.Ictal);
            
            for j = 1:length(pt(whichPt).sz)
                for k = 1:length(szs)
                    sz = info.Events.Ictal.(szs{k});
                    if isfield(sz,'SeizureEEC') == 0
                       continue
                    end
                    
                    if pt(whichPt).sz(j).onset == sz.SeizureEEC
                        
                        s = regexp(sz.FILE,'-\d.mat');
                        s = sz.FILE(s+1);
                        pt(whichPt).sz(j).adj_num = str2num(s);
                        
                        if isequal(pt(whichPt).sz(j).electrodes,...
                                sz.SEIZURE_ONSET_ELECTRODES) == 0
                            error('what\n');
                        end
                    end
                end
            end
        end
    end
end


end