function cc_comparison(stats,pt)

%% Look at contiguous and -5 seconds
contig_text = 'contiguous';
sec_text = 'sec_neg5';

std_diff_cc = [];
std_diff_cc_regional = [];

for i = 1:length(stats)
    
    locs = pt(i).new_elecs.locs;
    
    %% Extract just numbers from name (for plotting)
    %names = [names;stats(i).name];
    %[num_idx_s] = regexp(stats(i).name,'\d+');
    %name_nums = [name_nums;stats(i).name(num_idx_s:end)];
    
    % What is the min_cc
    temp_min_cc = stats(i).cc.(contig_text).(sec_text).min.min_elec;
    %min_cc = [min_cc;temp_min_cc];
    
    % What is the min region
    temp_min_elecs = stats(i).cc.(contig_text).(sec_text).regional_cc.min_elecs;
    
    
    
    % Plot
    if 1 == 0
    figure
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
    hold on
    if isnan(temp_min_elecs) == 0
        scatter3(locs(temp_min_elecs,1),locs(temp_min_elecs,2),locs(temp_min_elecs,3),100,'b');
    end
    scatter3(locs(temp_min_cc,1),locs(temp_min_cc,2),locs(temp_min_cc,3),60,'r','filled');
    pause
    close(gcf)
    end
    
    
    % What is the standard deviation of the distance between true min cc
    % and min cc in resampled network
    
    temp_std_diff_cc = stats(i).cc.(contig_text).(sec_text).min.dist_std(4);
    std_diff_cc = [std_diff_cc;temp_std_diff_cc];
    
    % What is the standard deviation of the distance between true min cc
    % regional and min cc regional in resampled network
    temp_std_diff_cc_regional = stats(i).cc.(contig_text).(sec_text).regional_cc.dist_std(4);
    std_diff_cc_regional = [std_diff_cc;temp_std_diff_cc_regional];
    
end

end