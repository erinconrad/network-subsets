function cc_comparison(stats,pt)

%% Look at contiguous and -5 seconds
contig_text = 'random';
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
    temp_min_cc = stats(i).cc.(contig_text).(sec_text).min_cc.true;
    %min_cc = [min_cc;temp_min_cc];
    
    % What is the min region
    temp_min_elecs = stats(i).cc.(contig_text).(sec_text).regional_cc.true;
    if isnan(temp_min_elecs) == 1
        centroid = nan;
    else
        centroid = mean(locs(temp_min_elecs,:),1);
    end
    
    % What is the standard deviation of the min cc loc in the resampled
    % network
    
    temp_std_diff_cc = stats(i).cc.(contig_text).(sec_text).min_cc.res_std(4,:);
    std_diff_cc = [std_diff_cc;temp_std_diff_cc];
    
    % What is the standard deviation of min cc region centroid in the
    % resampled network
    temp_std_diff_cc_regional = stats(i).cc.(contig_text).(sec_text).regional_cc.res_std(4,:);
    std_diff_cc_regional = [std_diff_cc;temp_std_diff_cc_regional];
    
    
    % Plot
    if 1 == 1
    n = 30;
    figure
    set(gcf,'Position',[200 200 1000 400])
    subplot(1,2,1)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
    hold on
    if isnan(temp_min_elecs) == 0
        scatter3(locs(temp_min_elecs,1),locs(temp_min_elecs,2),locs(temp_min_elecs,3),...
            100,'b','filled');
    end
    scatter3(locs(temp_min_cc,1),locs(temp_min_cc,2),locs(temp_min_cc,3),60,'r','filled');
    xticklabels([])
    yticklabels([])
    zticklabels([])
    title(sprintf(['Location of min control centrality\nelectrode '...
        'and region for %s'],pt(i).name));
    set(gca,'fontsize',20)
    
    subplot(1,2,2)
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
    hold on
    scatter3(locs(temp_min_cc,1),locs(temp_min_cc,2),locs(temp_min_cc,3),60,'r','filled');
    [x,y,z] = ellipsoid(locs(temp_min_cc,1),locs(temp_min_cc,2),locs(temp_min_cc,3),...
        temp_std_diff_cc(1),temp_std_diff_cc(2),temp_std_diff_cc(3),n);
    C(:,:,1) = ones(n); C(:,:,2) = zeros(n); C(:,:,3) = zeros(n);
    p_s = surf(x,y,z,C,'EdgeColor','none');
    alpha(p_s,0.2);
    if isnan(temp_min_elecs) == 0
        scatter3(centroid(1),centroid(2),centroid(3),100,'b','filled');
        [x,y,z] = ellipsoid(centroid(1),centroid(2),centroid(3),temp_std_diff_cc_regional(1),...
        temp_std_diff_cc_regional(2),temp_std_diff_cc_regional(3),n);
        C(:,:,1) = zeros(n); C(:,:,2) = zeros(n); C(:,:,3) = ones(n);
        p_r = surf(x,y,z,C,'EdgeColor','none');
        alpha(p_r,0.2);
        if isempty(pt(i).resec) == 0
            resec = pt(i).resec.nums;
            scatter3(locs(resec,1),locs(resec,2),locs(resec,3),'*')
        end
    end
    xticklabels([])
    yticklabels([])
    zticklabels([])
    title(sprintf(['Location and std of min control centrality\nelectrode '...
        'and region centroid for %s'],pt(i).name));    
    set(gca,'fontsize',20)
    
   % fprintf('Std of min cc loc is %1.1f mm\n',std_diff_cc);
   % fprintf('Std of min cc region centroid is %1.1f mm\n',std_diff_cc);
    pause
    close(gcf)
    end
    
    
    
    
end

end