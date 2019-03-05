function new_cc_comparison(stats,pt)

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'basic_metrics/'];

%% Look at contiguous and -5 seconds
contig_text = 'random';
sec_text = 'sec_neg5';

np = length(stats);

%% Get data
single_true = zeros(np,1);
region_true =cell(np,1);
single_95 = cell(np,1);
region_95 = cell(np,1);
single_names = cell(np,1);
region_names = cell(np,1);
resec = cell(np,1);
ds = zeros(np,1);
good_outcome = zeros(np,1);


for i = 1:length(stats)
    single_true(i) = stats(i).(contig_text).(sec_text).min_cc.true';
    single_names{i} = pt(i).name;
    single_95{i} = stats(i).(contig_text).(sec_text).min_cc_elecs.single_95';
    nchs = length(pt(i).new_elecs.electrodes);
    
    outcome = pt(i).clinical.outcome;
    if contains(outcome,'1.') == 1 || contains(outcome,'ILAE1') == 1||...
            contains(outcome,'ILAE2') == 1
        good_outcome_temp = 1;
    else
        good_outcome_temp = 0;
    end
    good_outcome(i) = good_outcome_temp;        
    
    
    if sum(isnan(stats(i).(contig_text).(sec_text).regional_cc.true)) == 0
        resec_temp = pt(i).resec.nums;
        resec{i} = resec_temp;
        region_true{i} = stats(i).(contig_text).(sec_text).regional_cc.true';
        region_names{i} = pt(i).name;
        region_95{i} = stats(i).(contig_text).(sec_text).min_cc_elecs.regional_95';
        is_resec = zeros(nchs,1); is_resec(resec_temp) = 1;
        is_sync = zeros(nchs,1); is_sync(stats(i).(contig_text).(sec_text).regional_cc.true') = 1;
        
        ds_temp = dice_score(is_resec,is_sync);
        ds(i) = ds_temp;
    else
        ds(i) = nan;
        resec{i} = nan;
        region_true{i} = nan;
        region_names{i} = nan;
        region_95{i} = nan;
        
    end

end

%% Do stats

% Is the dice score for good outcome better than for bad outcome?
[p,h,stats] = ranksum(ds(good_outcome==1),ds(good_outcome==0));
fprintf(['Mean dice score of good outcome is %1.2f and bad outcome %1.2f.\n'...
    'Wilcoxon rank sum: p = %1.2f, ranksum = %d.\n'],nanmean(ds(good_outcome==1)),...
    nanmean(ds(good_outcome==0)),p,stats.ranksum);

% Get the number of electrodes in the 95% CI for single min cc
num_single_95 = cellfun(@length,single_95);

fprintf(['The median number of electrodes in the 95%% CI for the most\n'...
    'synchronizing electrode is %1.1f, range %d-%d\n'],nanmedian(num_single_95),...
    min(num_single_95),max(num_single_95));

% Get the number of electrodes in the 95% CI for regional min cc
num_region_95 = cellfun(@(x) sum(~isnan(unique((x)))),region_95);
num_true_region = cellfun(@(x) sum(~isnan(unique((x)))),region_true);
mult_region_95 = num_region_95./num_true_region;

fprintf(['The mean number of electrodes in the 95%% CI for the most\n'...
    'synchronizing region is %1.1f, range %d-%d\n'],nanmean(num_region_95),...
    min(num_region_95),max(num_region_95));

fprintf(['The mean ratio of electrodes in the 95%% CI for the most\n'...
    'synchronizing region to true synchronizing region is %1.1f, range %1.1f-%1.1f\n'],nanmean(mult_region_95),...
    min(mult_region_95),max(mult_region_95));


%% Make plots
%% Individual plots
for i = 1:length(stats)
    locs = pt(i).new_elecs.locs;
    
    if isempty(pt(i).resec) == 1, continue; end
    resec = locs(pt(i).resec.nums,:);
    name = pt(i).name;
    pt_folder = [outFolder,name,'/'];
    ex_regional = [];
    true = stats(i).(contig_text).(sec_text).regional_cc.true';
    ex_regional = [ex_regional;true];
    
    
    figure
    c = colormap(parula(5));
    pl = zeros(6,1);
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
    hold on
    count = 1;
    col_one = ones(size(true,1),1).*c(1,:);
    pl(1)=scatter3(locs(true,1),locs(true,2),locs(true,3),...
            100,col_one,'filled');
    
    for j = [70 80 90 95]
        count = count + 1;
        elecs = stats(i).(contig_text).(sec_text).min_cc_elecs.(sprintf('regional_%d',j))';
        elecs(ismember(elecs,ex_regional)) = [];
        ex_regional = [ex_regional;elecs];
        c=colormap(parula(5));
        col_all = ones(size(elecs,1),1).*c(count,:);  
        pl(count)=scatter3(locs(elecs,1),locs(elecs,2),locs(elecs,3),...
                100,col_all,'filled');
    end
    
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
        pl(6) = scatter3(resec(:,1),resec(:,2),resec(:,3),15,'k','filled');
    legend(pl,{'True','70%','80%','90%','95%','Resected'},'Position',[0.88 0.7 0.03 0.2],...
                'box','off');
    title(sprintf('%s, dice score = %1.2f',name,ds(i)))
    set(gca,'fontsize',20)

    xticklabels([])
    yticklabels([])
    zticklabels([])
    pause
    print(gcf,[pt_folder,'most_sync_',contig_text,sec_text],'-depsc');
    close(gcf)
end

%% 3 patient plot
figure
set(gcf,'Position',[100 100 1100 600]);
[ha,pos] = tight_subplot(2, 3, [0.01 0.03], [0.02 0.07],[0.12 0.08]);
%delete(ha(3))
count = 0;

for text = {'single','regional'}
    for i = [1 4 7]
        count = count + 1;
        %if count == 3, count =  count+ 1; end

        locs = pt(i).new_elecs.locs;
        resec = locs(pt(i).resec.nums,:);
        axes(ha(count))
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
        hold on
        new_count = 0;
        pl = zeros(6,1);
        c = colormap(parula(5));
        ex_single = [];
        ex_regional = [];
        if strcmp(text,'single') == 1
            true1 = stats(i).(contig_text).(sec_text).min_cc.true';
            ex_single = [ex_single;true1];
        else
            true1 = stats(i).(contig_text).(sec_text).regional_cc.true';
            ex_regional = [ex_regional;true1];
        end
        col_one = ones(size(true1,1),1).*c(1,:);
        pl(1)=scatter3(locs(true1,1),locs(true1,2),locs(true1,3),...
                100,col_one,'filled');
        
        for j = [70 80 90 95]
            new_count = new_count + 1;
            elecs = stats(i).(contig_text).(sec_text).min_cc_elecs.([text{1},sprintf('_%d',j)])';
            if strcmp(text,'single') == 1
                elecs(ismember(elecs,ex_single)) = [];
                ex_single = [ex_single;elecs];
                c=colormap(parula(5));
                col_all = ones(size(elecs,1),1).*c(new_count+1,:);
                
            else
                elecs(ismember(elecs,ex_regional)) = [];
                ex_regional = [ex_regional;elecs];
                c=colormap(parula(5));
                col_all = ones(size(elecs,1),1).*c(new_count+1,:);   
            end

            if sum(isnan(elecs)) > 0
                continue
            end


           pl(new_count+1)=scatter3(locs(elecs,1),locs(elecs,2),locs(elecs,3),...
                100,col_all,'filled');
           
 
        end
        
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
        pl(6) = scatter3(resec(:,1),resec(:,2),resec(:,3),15,'k','filled');
        if count == 3
            legend(pl,{'True','70%','80%','90%','95%','Resected'},'Position',[0.91 0.75 0.05 0.2],...
                'box','on');
        elseif count == 6
            legend(pl,{'True','70%','80%','90%','95%','Resected'},'Position',[0.91 0.27 0.05 0.2],...
                'box','on');
        end

        
        annotation('textbox',[0 0.7 0.1 0.1],'String',...
    sprintf('Most\nsynchronizing\nelectrode'),'LineStyle','none','fontsize',20);

        annotation('textbox',[0 0.22 0.1 0.1],'String',...
    sprintf('Most\nsynchronizing\nregion'),'LineStyle','none','fontsize',20);

        annotation('textbox',[0.18 0.9 0.1 0.1],'String',...
    sprintf('Patient 1'),'LineStyle','none','fontsize',20);

        annotation('textbox',[0.48 0.9 0.1 0.1],'String',...
    sprintf('Patient 2'),'LineStyle','none','fontsize',20);


        annotation('textbox',[0.77 0.9 0.1 0.1],'String',...
    sprintf('Patient 3'),'LineStyle','none','fontsize',20);
        
        set(gca,'fontsize',20)

        xticklabels([])
        yticklabels([])
        zticklabels([])
%}
    
    end

end

print(gcf,[outFolder,'pc_95_cc_',contig_text,sec_text],'-depsc');
print(gcf,[outFolder,'pc_95_cc_',contig_text,sec_text],'-dpng');



end