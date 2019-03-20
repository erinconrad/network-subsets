function socc = second_order_cc(A)

nchs = size(A,1);

% Calculate the control centrality
cc = control_centrality(A);

% calculate resampled control centralities, removing a single electrode at
% a time
cc_resample = nan(nchs,nchs);

for i = 1:nchs
    A_temp = A;
    
    % Remove a single electrode
    A_temp(i,:) = [];
    A_temp(:,i) = [];
    curr_chs = 1:nchs;
    curr_chs(i) = [];
    
    % Recalculate cc of remaining electrodes
    cc_temp = control_centrality(A_temp);
    
    % Fill it in with the appropriate electrodes
    for j = 1:length(curr_chs)
        cc_resample(i,curr_chs(j)) = cc_temp(j);
    end
    
    
end

%% Calculate agreement between resampled and original
agree = nan(nchs,1);
for i = 1:nchs
    new_cc = cc_resample(i,:);
    agree(i) = corr(new_cc',cc,'Type','Spearman','Rows','complete');
end

socc = agree;

end