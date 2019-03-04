function new_cc_comparison(stats,pt)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'basic_metrics/'];

%% Look at contiguous and -5 seconds
contig_text = 'contiguous';
sec_text = 'sec_neg5';

figure
set(gcf,'Position',[100 100 1000 600]);
[ha,pos] = tight_subplot(2, 3, [0.1 0.03], [0.02 0.14],[0.05 0.01]);
%delete(ha(3))
count = 0;


for i = [1 2 4 8 9 33]
    count = count + 1;
    %if count == 3, count =  count+ 1; end
    
    locs = pt(i).new_elecs.locs;
    temp_min_cc = stats(i).(contig_text).(sec_text).min_cc.true;
    temp_min_elecs = stats(i).(contig_text).(sec_text).regional_cc.true;
    elecs_95 = stats(i).(contig_text).(sec_text).min_cc.elecs_95';
    regional_95 = stats(i).(contig_text).(sec_text).regional_cc.elecs_95';
    if sum(isnan(regional_95)) > 0
        continue
    end

    axes(ha(count))
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
    hold on
    scatter3(locs(regional_95,1),locs(regional_95,2),locs(regional_95,3),...
        100,ones(size(regional_95,1),1).*[0.7 0.7 1],'filled');
%{
    scatter3(locs(elecs_95,1),locs(elecs_95,2),locs(elecs_95,3),...
        100,ones(size(elecs_95,1),1).*[1 0.7 0.7],'filled');
%}
    scatter3(locs(temp_min_elecs,1),locs(temp_min_elecs,2),locs(temp_min_elecs,3),...
        60,'b','filled');
%{
    scatter3(locs(temp_min_cc,1),locs(temp_min_cc,2),locs(temp_min_cc,3),...
        60,'r','filled');
%}
    

end


end