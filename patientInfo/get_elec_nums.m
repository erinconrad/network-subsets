function get_elec_nums(pt,stats)

names = {};
num_elecs = [];
num_soz = [];
patient_id = [];

for i = 1:33
    if isempty(stats(i).name) == 1, continue; end
    if strcmp(stats(i).name,pt(i).name) == 0, error('whatn'); end
    names = [names;pt(i).name];
    patient_id = [patient_id;i];
    num_elecs = [num_elecs;length(pt(i).new_elecs.electrodes)];
    num_soz = [num_soz;length(pt(i).soz.nums)];
    
end

num_non_soz = num_elecs - num_soz;
poss_combs = nan(size(num_soz));

for i = 1:length(num_soz)
    if num_non_soz(i) < num_soz(i)
        poss_combs(i) = nan;
    else
        poss_combs(i) = nchoosek(num_non_soz(i),num_soz(i));
    end
end

table(names,patient_id,num_elecs,num_non_soz,poss_combs)

%{ 
Concerning patients:
- HUP064 (1)
- HUP078 (8)
- HUP116 (21)
- Study012 (23)
- Study017 (25)
- Study022 (29)
%}

end