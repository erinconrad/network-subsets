function cc_comparison(stats,pt)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'basic_metrics/'];

%% Look at contiguous and -5 seconds
contig_text = 'random';
sec_text = 'sec_neg5';



names_cc = {};
names_cc_regional = {};
std_diff_cc_regional = [];
std_diff_cc = [];
main_cc_std_80 = [];
cc_95 = [];
cc_95_regional = [];

for i = 1:length(stats)
    temp_std_diff_cc = stats(i).(contig_text).(sec_text).min_cc.res_std(4,:);
    if isfield(stats(i).(contig_text).(sec_text).min_cc,'pct_95') == 0
        continue
    end
    temp_95_diff_cc = diff(squeeze(stats(i).(contig_text).(sec_text).min_cc.pct_95(4,:,:)),1)/2;
    if isequal(temp_std_diff_cc,[0 0 0]) == 0
        names_cc =[names_cc;pt(i).name];
        std_diff_cc = [std_diff_cc;temp_std_diff_cc];
        cc_95 = [cc_95;temp_95_diff_cc];
    end
    
    temp_std_diff_cc_regional = stats(i).(contig_text).(sec_text).regional_cc.res_std(4,:);
    temp_95_diff_cc_regional = diff(squeeze(stats(i).(contig_text).(sec_text).regional_cc.pct_95(4,:,:)),1)/2;
    if isequaln(temp_std_diff_cc_regional,[nan nan nan]) == 0
        names_cc_regional =[names_cc_regional;pt(i).name];
        std_diff_cc_regional = [std_diff_cc_regional;temp_std_diff_cc_regional];
        cc_95_regional = [cc_95_regional;temp_95_diff_cc_regional];
    end
    
    main_cc_std_80 = [main_cc_std_80;stats(i).(contig_text).(sec_text).cc.rel_std(4)];
end

std_diff_cc = vecnorm(std_diff_cc,2,2);
std_diff_cc_regional = vecnorm(std_diff_cc_regional,2,2);
cc_95 = vecnorm(cc_95,2,2);
cc_95_regional = vecnorm(cc_95_regional,2,2);

fprintf('Range of 95%% of min cc is %1.1f-%1.1f\n.',min(cc_95),max(cc_95));
fprintf('Range of 95%% of min cc region is %1.1f-%1.1f\n.',min(cc_95_regional),max(cc_95_regional));

table(names_cc,cc_95)
table(names_cc_regional,cc_95_regional)

figure
set(gcf,'Position',[100 100 1000 600]);
[ha,pos] = tight_subplot(2, 3, [0.1 0.03], [0.02 0.14],[0.05 0.01]);
%delete(ha(3))
count = 0;


for i = [1 4 8 32 33]
    count = count + 1;
    %if count == 3, count =  count+ 1; end
    
    locs = pt(i).new_elecs.locs;
    
    %% Extract just numbers from name (for plotting)
    %names = [names;stats(i).name];
    %[num_idx_s] = regexp(stats(i).name,'\d+');
    %name_nums = [name_nums;stats(i).name(num_idx_s:end)];
    
    % What is the min_cc
    temp_min_cc = stats(i).(contig_text).(sec_text).min_cc.true;
    %min_cc = [min_cc;temp_min_cc];
    
    % What is the min region
    temp_min_elecs = stats(i).(contig_text).(sec_text).regional_cc.true;
    if isnan(temp_min_elecs) == 1
        centroid = nan;
    else
        centroid = mean(locs(temp_min_elecs,:),1);
    end
    
    % What is the standard deviation of the min cc loc in the resampled
    % network
    
    temp_std_diff_cc = stats(i).(contig_text).(sec_text).min_cc.res_std(4,:);
    %std_diff_cc = [std_diff_cc;temp_std_diff_cc];
    
    % What are the half-lengths of the 95th CIs of the min cc loc in the
    % resampled network?
    temp_95_diff_cc = diff(squeeze(stats(i).(contig_text).(sec_text).min_cc.pct_95(4,:,:)),1)/2;
    single_middle = mean(squeeze(stats(i).(contig_text).(sec_text).min_cc.pct_95(4,:,:)),1);
    
    % What is the standard deviation of min cc region centroid in the
    % resampled network
    temp_std_diff_cc_regional = stats(i).(contig_text).(sec_text).regional_cc.res_std(4,:);
    %std_diff_cc_regional = [std_diff_cc_regional;temp_std_diff_cc_regional];
    
    temp_95_diff_cc_regional = diff(squeeze(stats(i).(contig_text).(sec_text).regional_cc.pct_95(4,:,:)),1)/2;
    regional_middle = mean(squeeze(stats(i).(contig_text).(sec_text).regional_cc.pct_95(4,:,:)),1);
    
    
    n = 30;
    % Plot
    if i == 1
    
        axes(ha(count))
        count = count + 1;
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
        title(sprintf(['Location of most synchronizing\nelectrode '...
            'and region for %s'],pt(i).name));
        set(gca,'fontsize',20)
    end
    
    axes(ha(count))
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
    hold on
    scatter3(locs(temp_min_cc,1),locs(temp_min_cc,2),locs(temp_min_cc,3),60,'r','filled');
    
    [x,y,z] = ellipsoid(locs(temp_min_cc,1),locs(temp_min_cc,2),locs(temp_min_cc,3),...
        temp_95_diff_cc(1),temp_95_diff_cc(2),temp_95_diff_cc(3),n);
    %{
    [x,y,z] = ellipsoid(single_middle(1),single_middle(2),single_middle(3),...
        temp_95_diff_cc(1),temp_95_diff_cc(2),temp_95_diff_cc(3),n);
        %}
    C(:,:,1) = ones(n); C(:,:,2) = zeros(n); C(:,:,3) = zeros(n);
    p_s = surf(x,y,z,C,'EdgeColor','none');
    alpha(p_s,0.2);
    if isnan(temp_min_elecs) == 0
        scatter3(centroid(1),centroid(2),centroid(3),100,'b','filled');
        
        [x,y,z] = ellipsoid(centroid(1),centroid(2),centroid(3),temp_95_diff_cc_regional(1),...
        temp_95_diff_cc_regional(2),temp_95_diff_cc_regional(3),n);
        %{
        [x,y,z] = ellipsoid(regional_middle(1),regional_middle(2),regional_middle(3),temp_95_diff_cc_regional(1),...
        temp_95_diff_cc_regional(2),temp_95_diff_cc_regional(3),n);
        %}
        C(:,:,1) = zeros(n); C(:,:,2) = zeros(n); C(:,:,3) = ones(n);
        p_r = surf(x,y,z,C,'EdgeColor','none');
        alpha(p_r,0.2);
        %{
        if isempty(pt(i).resec) == 0
            resec = pt(i).resec.nums;
            scatter3(locs(resec,1),locs(resec,2),locs(resec,3),'*')
        end
        %}
    end
    xticklabels([])
    yticklabels([])
    zticklabels([])
    if count == 2
        title(sprintf(['95%% confidence intervals\n'...
            'for %s'],pt(i).name));
    else
        title(sprintf('%s',pt(i).name));
    end
    set(gca,'fontsize',20)
    print(gcf,[outFolder,'pc_95_cc_',contig_text,sec_text],'-depsc');
    
   % fprintf('Std of min cc loc is %1.1f mm\n',std_diff_cc);
   % fprintf('Std of min cc region centroid is %1.1f mm\n',std_diff_cc);
   
    
end

end