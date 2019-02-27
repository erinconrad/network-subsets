function [names,nums] = getElecNums(pt,whichPt,names)
nums = [];
discard = [];

for i = 1:length(names)
    foundit = 0;
    for j = 1:length(pt(whichPt).new_elecs.electrodes)
        if strcmp(names{i},pt(whichPt).new_elecs.electrodes(j).name) == 1
            nums = [nums;j];
            foundit = 1;
            break
        end
    end
    if foundit == 0
        fprintf('Warning, did not find %s for %s\n',names{i},pt(whichPt).name);
        discard = [discard;i];
    end
end

names(discard) = [];

end