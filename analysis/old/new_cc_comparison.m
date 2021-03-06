function new_cc_comparison(stats,pt)

%% Parameters
contig_text = 'random';
sec_text = 'sec_0';
freq = 'high_gamma';
do_individ_plots = 0;
which_texts = {'ns','min_cc_elecs'};
whichPts = [1 10 21];

%% Locations

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'cc_comparison/'];

np = length(stats);

%% Initialize arrays
single_true = nan(np,1);
region_true =cell(np,1);
single_95 = cell(np,1);
region_95 = cell(np,1);
single_names = cell(np,1);
region_names = cell(np,1);
resec = cell(np,1);
ds = nan(np,1);
good_outcome = nan(np,1);
min_reg_cc = cell(np,1);

% Loop through patients
for i = 1:length(stats)
    
    if isempty(stats(i).name) ==1 ,continue; end
    
    base = stats(i).(freq).(contig_text).(sec_text);
    
    % Get the most synchronizing electrode
    single_true(i) = base.min_cc.true';
    single_names{i} = pt(i).name;
    single_95{i} = base.min_cc_elecs.single_95';
    nchs = length(pt(i).new_elecs.electrodes);
    
    % Get the outcome
    outcome = pt(i).clinical.outcome;
    if contains(outcome,'1.') == 1 || contains(outcome,'ILAE1') == 1||...
            contains(outcome,'ILAE2') == 1
        good_outcome_temp = 1;
    else
        good_outcome_temp = 0;
    end
    good_outcome(i) = good_outcome_temp;        
    
    % Only do analysis if there were resected electrodes
    if sum(isnan(base.regional_cc.true)) == 0
        
        % Get resected electrodes
        resec_temp = pt(i).resec.nums;
        resec{i} = resec_temp;
        n_res = length(resec_temp);
        
        % Get the electrodes in the most synchronizing region
        region_true{i} = base.regional_cc.true';
        region_names{i} = pt(i).name;
        region_95{i} = base.min_cc_elecs.regional_95';
        
        % Get electrodes with minimal regional control centrality
        %{
        reg_cc = base.cc_reg.true;
        [~,reg_cc_sorted_chs] = sort(reg_cc);
        min_reg_cc{i} = reg_cc_sorted_chs(1:n_res);
        %}
        
        % Get dice score to examine overlap between resected electrodes and
        % most synchronizing electrodes
        is_resec = zeros(nchs,1); is_resec(resec_temp) = 1;
        is_sync = zeros(nchs,1); is_sync(region_true{i}') = 1;
        if sum(is_resec) ~= sum(is_sync), error('what\n'); end
        ds_temp = dice_score(is_resec,is_sync);
        ds(i) = ds_temp;
        
    else
        ds(i) = nan;
        resec{i} = nan;
        region_true{i} = nan;
        region_names{i} = nan;
        region_95{i} = nan;
        min_reg_cc{i} = nan;
        
    end

end

%% Do stats
if strcmp(which_texts{1},'ns') == 1
    box_text = 'node strength';
else
    box_text = '???';
end

%{
% Is the dice score for good outcome better than for bad outcome?
[p,~,stats1] = ranksum(ds(good_outcome==1),ds(good_outcome==0));
fprintf(['Mean dice score of good outcome is %1.2f and bad outcome %1.2f.\n'...
    'Wilcoxon rank sum: p = %1.2f, ranksum = %d.\n'],nanmean(ds(good_outcome==1)),...
    nanmean(ds(good_outcome==0)),p,stats1.ranksum);

% Do a chi2 comparing good and bad outcome with overlap and no overlap
yes_overlap = (ds>0);
no_overlap = (ds==0);
yes_overlap_outcome = good_outcome(yes_overlap);
no_overlap_outcome = good_outcome(no_overlap);

[tbl,chi2,p,labels] = crosstab([ones(length(yes_overlap_outcome),1);...
                2*ones(length(no_overlap_outcome),1)],...
                [yes_overlap_outcome;no_overlap_outcome]);

fprintf(['Chi2 comparing number of people with good outcome among those\n'...
    'with overlap and those without: p = %1.3f, chi2 = %1.2f\n'],p,chi2);

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

%% Compare number of electrodes and ratio of electrodes in 95% CI with outcome

% Number
[p,~,stats1] = ranksum(num_region_95(good_outcome==1),num_region_95(good_outcome==0));
fprintf(['Mean number of electrodes in 95%% CI for good outcome is %1.2f and bad outcome %1.2f.\n'...
    'Wilcoxon rank sum: p = %1.2f, ranksum = %d.\n'],nanmean(num_region_95(good_outcome==1)),...
    nanmean(num_region_95(good_outcome==0)),p,stats1.ranksum);

% Ratio
[p,~,stats1] = ranksum(mult_region_95(good_outcome==1),mult_region_95(good_outcome==0));
fprintf(['Mean ratio of electrodes in 95%% CI for good outcome is %1.2f and bad outcome %1.2f.\n'...
    'Wilcoxon rank sum: p = %1.2f, ranksum = %d.\n'],nanmean(mult_region_95(good_outcome==1)),...
    nanmean(mult_region_95(good_outcome==0)),p,stats1.ranksum);
%}

%% Make plots
%% Individual plots
if do_individ_plots == 1
for i = 1:length(stats)
    base = stats(i).(freq).(contig_text).(sec_text);
    locs = pt(i).new_elecs.locs;
    
    if isempty(pt(i).resec) == 1, continue; end
    resec = locs(pt(i).resec.nums,:);
    name = pt(i).name;
    pt_folder = [outFolder,name,'/'];
    ex_regional = [];
    true = base.regional_cc.true';
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
        elecs = base.min_cc_elecs.(sprintf('regional_%d',j))';
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
    %print(gcf,[pt_folder,'most_sync_',contig_text,sec_text],'-depsc');
    close(gcf)
end
end

%% 3 patient plot
figure
set(gcf,'Position',[100 100 1100 700]);
[ha,pos] = tight_subplot(3, 3, [0.01 0.03], [0.02 0.07],[0.12 0.08]);
del = 0.02;
%{
set(ha(7),'position',[pos{4}(1), pos{7}(2), pos{5}(1) - pos{4}(1) + ...
    pos{5}(3)/2 - del, pos{7}(4)]);
set(ha(8),'position',[pos{5}(1) + pos{5}(3)/2 + del, pos{7}(2),...
    pos{5}(1) - pos{4}(1) + pos{5}(3)/2 - del, pos{7}(4)]);
%}
set(ha(7),'position',[pos{4}(1), pos{7}(2), (pos{5}(1) - pos{4}(1) + ...
    pos{5}(3)/2)*2, pos{7}(4)]);
delete(ha(8))

delete(ha(9))
%delete(ha(3))
count = 0;
text_count = 0;

for text = which_texts
    text_count = text_count + 1;
    for i = whichPts
        base = stats(i).(freq).(contig_text).(sec_text);
        count = count + 1;
        %if count == 3, count =  count+ 1; end

        locs = pt(i).new_elecs.locs;
        resec = locs(pt(i).resec.nums,:);
        axes(ha(count))
        
        %hold on
        new_count = 0;
        pl = zeros(6,1);
        c = colormap(parula(5));
        ex_single = [];
        ex_regional = [];
        if strcmp(text,'min_cc_elecs') == 0
            true1 = base.(text{1}).true';
            ex_single = [ex_single;true1];
        else
            true1 = base.regional_cc.true';
            ex_regional = [ex_regional;true1];
        end
        col_one = ones(size(true1,1),1).*c(1,:);
        pl(1)=scatter3(locs(true1,1),locs(true1,2),locs(true1,3),...
                100,col_one,'filled');
        hold on
        
        for j = [70 80 90 95]
            new_count = new_count + 1;
            
            if strcmp(text,'min_cc_elecs') == 0
                elecs = base.(text{1}).(sprintf('single_%d',j))';
                elecs(ismember(elecs,ex_single)) = [];
                ex_single = [ex_single;elecs];
                c=colormap(parula(5));
                col_all = ones(size(elecs,1),1).*c(new_count+1,:);
                
            else
                elecs = base.min_cc_elecs.(sprintf('regional_%d',j))';
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
        
        if text_count == 1, title(sprintf('%s',stats(i).name)); end
        
        pl(6) = scatter3(resec(:,1),resec(:,2),resec(:,3),15,'k','filled');
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
        if count == 3
            
        elseif count == 6
            legend(pl,{'True','70%','80%','90%','95%','Resected'},'Position',[0.91 0.75 0.05 0.2],...
                'box','on');
        end

        if count == 1 || count == 4
            view(150,32.4);
        elseif count == 2 || count == 5
            view(-94,14)
        elseif count == 3 || count == 6
            view(-31.1,6.8)
        end
       
        set(gca,'fontsize',20)
        xticklabels([])
        yticklabels([])
        zticklabels([])
        
        
        if count == 1 
            zlabel(sprintf('Highest\n%s',box_text));
        elseif count == 4
            zlabel({'Highest regional','control centrality'});
        end
        %}
%}
    
    end
    %{
    %% Make the histograms
    if text_count == 1
        axes(ha(7))
    else
        axes(ha(8))
    end
    rat = [];
    
    for i = 1:length(stats)
        if isfield(stats(i),(freq)) == 0, continue; end
        if isempty(stats(i).(freq)) == 1, continue; end
        base = stats(i).(freq).(contig_text).(sec_text); 
        if strcmp(text,'min_cc_elecs') == 0
            rat_t = length(base.(text{1}).single_95);
            rat = [rat;rat_t];
        else
            if isempty(pt(i).resec) == 1, rat = [rat;nan]; continue; end
            rat_t = length(base.min_cc_elecs.regional_95);%/length(base.regional_cc.true);
            rat= [rat;rat_t];
        end
        
    end
    
    bar(rat)
    set(gca,'ylim',[0 max(rat)]);
    xticklabels([])
    xlabel('Which patient');
    set(gca,'fontsize',20);
    %}
    
end
    
axes(ha(7))
all_rat = cell(6,1);
all_mets = {'ns','bc','ec','clust','cc'};
all_met_names = {'Node strength','Betweenness centrality',...
    'Eigenvector centrality','Clustering coefficient',...
    'Control centrality','Regional control centrality'};
for i = 1:length(stats)
    if isempty(stats(i).(freq)) == 1, continue; end
    for metric = 1:length(all_mets)
    
        base = stats(i).(freq).(contig_text).(sec_text);
        
        all_rat{metric} = [all_rat{metric};...
            length(base.(all_mets{metric}).single_95)];
        
    end
    
    if isempty(pt(i).resec) == 1, continue; end
    all_rat{6} = [all_rat{6};...
        length(base.min_cc_elecs.regional_95)/...
        length(base.regional_cc.true)];
end

for i = 1:length(all_rat)
    scatter(i*ones(size(all_rat{i},1),1)+0.05*randn(size(all_rat{i},1),1)...
        ,all_rat{i},...
       80,'filled','MarkerEdgeColor',[0 0.4470 0.7410],...
       'MarkerFaceColor',[0 0.4470 0.7410])
    hold on
end
xticks(1:length(all_rat))
xticklabels(all_met_names)
fix_xticklabels(gca,0.1,{'FontSize',20});
ylabel({'Relative number of','electrodes in 95% CI','(all patients)'})



%{
 annotation('textbox',[0 0.75 0.1 0.1],'String',...
sprintf('Highest\n%s\nelectrode',box_text),'LineStyle','none','fontsize',20);

    annotation('textbox',[0 0.45 0.1 0.1],'String',...
sprintf('Most\nsynchronizing\nregion'),'LineStyle','none','fontsize',20);
%}

%{
    annotation('textbox',[0.18 0.9 0.1 0.1],'String',...
sprintf('Patient 1'),'LineStyle','none','fontsize',20);

    annotation('textbox',[0.48 0.9 0.1 0.1],'String',...
sprintf('Patient 2'),'LineStyle','none','fontsize',20);


    annotation('textbox',[0.77 0.9 0.1 0.1],'String',...
sprintf('Patient 3'),'LineStyle','none','fontsize',20);
%}

    set(gca,'fontsize',20)

print(gcf,[outFolder,'nodal_',freq,contig_text,sec_text],'-depsc');
%print(gcf,[outFolder,'nodal_cc_',freq,contig_text,sec_text],'-dpng');



end