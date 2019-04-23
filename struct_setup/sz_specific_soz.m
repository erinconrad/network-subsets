function pt = sz_specific_soz(pt)


for i = 1:length(pt)
    
    if isempty(pt(i).new_elecs) == 1
        continue
    end
    
    for j = 1:length(pt(i).sz)
        %% Get names
        soz_names = {};

        
        for k = 1:length(pt(i).sz(j).electrodes)
            soz_names = [soz_names,char(pt(i).sz(j).electrodes{k})];
        end
        

        soz_names = unique(soz_names)';

        %% Get nums
        [soz_names,soz_nums] = getElecNums(pt,i,soz_names);
        pt(i).sz(j).nums = soz_nums;
        pt(i).sz(j).names = soz_names;

        %% Double check that they line up
        for k = 1:length(pt(i).sz(j).nums)
            which_num = pt(i).sz(j).nums(k);
            if strcmp(pt(i).new_elecs.electrodes(which_num).name,pt(i).sz(j).names{k}) == 0
                error('Electrode names and nums do not line up for %s\n',...
                    pt(i).name);
            end
        end
    end
    
end


end