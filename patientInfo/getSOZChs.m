function pt = getSOZChs(pt)


for i = 1:length(pt)
    
    if isempty(pt(i).new_elecs) == 1
        continue
    end
    
    %% Get names
    soz_names = {};
    
    for j = 1:length(pt(i).sz)
        for k = 1:length(pt(i).sz(j).electrodes)
            soz_names = [soz_names,char(pt(i).sz(j).electrodes{k})];
        end
    end
    
    soz_names = unique(soz_names)';
    
    %% Get nums
    [soz_names,soz_nums] = getElecNums(pt,i,soz_names);
    
    pt(i).soz.names = soz_names;
    pt(i).soz.nums = soz_nums;
    
    %% Double check that they line up
    for j = 1:length(pt(i).soz.nums)
        which_num = pt(i).soz.nums(j);
        if strcmp(pt(i).new_elecs.electrodes(which_num).name,pt(i).soz.names{j}) == 0
            error('Electrode names and nums do not line up for %s\n',...
                pt(i).name);
        end
    end
    
end


end