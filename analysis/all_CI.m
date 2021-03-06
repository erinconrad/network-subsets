function all_CI(stats,pt,example)

%{
This function takes info about the number and identities of electrodes
forming the nth% CI for highest metric value and plots summary stuff about
them
%}

%% Parameters
contig_text = 'random';
sec_text = 'sec_0';
freq = 'high_gamma';
which_texts = {'ns','min_cc_elecs'};
which_names = {'node strength','regional control centrality'};
if example == 1
    whichPts = 1;
else
    whichPts = [1 10 21];
end
which_global = 'sync';
global_name = 'Synchronizability';
global_names_all = {'Synchronizability','Global efficiency','Transitivity'};

%% Locations
if example ~= 1
    [electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
    pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
    outFolder = [resultsFolder,'cc_comparison/'];
end

figure  
set(gcf,'Position',[100 0 900 1100]);
if example == 1
    for i = 1:12
        ha(i) = subplot(4, 3,i);
    end
    set(ha(1),'position',[0.1100    0.7750    0.2700    0.1850]);
    set(ha(4),'position',[0.1100    0.5200    0.2700    0.1850]);
    set(ha(7),'position',[0.1100    0.2650    0.4    0.1850]);
    set(ha(10),'position',[0.1200    0.0100    0.4250    0.1850]);
    set(ha(11),'position',[0.6650    0.0100    0.3150    0.1850]);
    delete(ha(2))
    delete(ha(3))
    delete(ha(5))
    delete(ha(6))
    delete(ha(8))
    delete(ha(9))
    delete(ha(12))
    
else
    [ha,pos] = tight_subplot(4, 3, [0.07 0.03], [0.01 0.04],[0.11 0.02]);
    set(ha(7),'position',[pos{4}(1), pos{7}(2), (pos{5}(1) - pos{4}(1) + ...
        pos{5}(3)/2)*2, pos{7}(4)]);
    delete(ha(8))
    delete(ha(9))
    set(ha(10),'position',[pos{4}(1)+0.01, pos{10}(2), (pos{5}(1) - pos{4}(1) + ...
        pos{5}(3)/2-0.01), pos{10}(4)]);
    set(ha(11),'position',[pos{4}(1) + (pos{5}(1) - pos{4}(1) + ...
        pos{5}(3)/2) + 0.12, pos{10}(2), (pos{5}(1) - pos{4}(1) + ...
        pos{5}(3)/2) - 0.12, pos{10}(4)]);
    delete(ha(12))
end

count = 0;
text_count = 0;

%% First three plots (first row) will be ns for 3 patients
% Second three plots (2nd row) will be regional cc for 3 patients

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
        pl = zeros(5,1);
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
        
       % pl(6) = scatter3(resec(:,1),resec(:,2),resec(:,3),15,'k','filled');
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
        

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
            zlabel(sprintf('Highest\n%s',which_names{text_count}));
        elseif count == 4
            zlabel({'Lowest regional','control centrality'});
        end
        %}
%}
    
    if example == 1, count = count + 2; end
    if count == 3
            
        elseif count == 6
            if example == 1
                l1 = legend(pl,{'True','70%','80%','90%','95%'},'Position',[0.37 0.72 0.03 0.10],...
                'box','on');
            else
           l1 = legend(pl,{'True','70%','80%','90%','95%'},'Position',[0.92 0.72 0.03 0.10],...
                'box','on');
            end
            pause(2) % do not delete
            for k = 1:length(l1.EntryContainer.NodeChildren)
                l1.EntryContainer.NodeChildren(k).Icon.Transform.Children.Children.Size = 14;
            end
        end
    
    end
    
end

%% The 7th plot (3rd row) is violin plot of sync distributions for all patients
axes(ha(7))
all_rho = zeros(1000,3);
new_count = 0;
all_95_ci_width = nan(29,3);
for j = {'sync','eff','trans'}
    new_count = new_count + 1;
    true = nan(29,1);
    all_perm = nan(29,1000);
    
    for i = 1:length(stats)    
        if isempty(stats(i).name) ==1
            true(i) = nan;
            all_perm(i,:) = nan(1,1000);
            continue; 
        end

        base = stats(i).(freq).(contig_text).(sec_text);

        true(i) = base.(j{1}).true;
        all_perm(i,:) = base.(j{1}).all(4,:);
    end
    inan = find(isnan(true));
    true(inan) = [];
    all_perm(inan,:) = [];
    for i = 1:length(all_rho)
        all_rho(i,new_count) = corr(true,all_perm(:,i),'Type','Spearman');
    end
    for i = 1:size(all_perm,1)
        w = prctile(all_perm(i,:),[2.5 97.5]);
        all_95_ci_width(i,new_count) = w(2)-w(1);
    end
    if strcmp(j{1},which_global) == 1
        violin(all_perm',true,'facecolor',[0 0.4470 0.7410])
        xticklabels([])
        ylabel(sprintf('%s',global_name));
        if example == 1
            
        else
            title('All patients');
        end
        set(gca,'fontsize',20)
    end
end

%% Next plot (4th row is summary stats for nodal metrics)
axes(ha(10))

%{
fc= [0 0.4470 0.7410;
    0.8500 0.3250 0.0980;
    0.9290 0.6940 0.1250;
    0.4940 0.1840 0.5560;
    0.4660 0.6740 0.1880;
    0.3010 0.7450 0.9330];
%}

fc = [0 0.4470 0.7410;0.8500 0.3250 0.0980;0.9290 0.6940 0.1250;...
        0.4940 0.1840 0.5560;0.4660 0.6740 0.1880;0.3010 0.7450 0.9330;...
        0.6350 0.0780 0.1840;0.75 0.5 0.5];

all_rat = cell(6,1);
all_mets = {'cc','ns','bc','ec','clust'};
all_met_names = {'Control centrality','Node strength','Betweenness centrality',...
    'Eigenvector centrality','Clustering coefficient',...
    'Regional control centrality'};
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
pl = zeros(length(all_rat),1);
for i = 1:length(all_rat)
    pl(i) = scatter(i*ones(size(all_rat{i},1),1)+0.05*randn(size(all_rat{i},1),1)...
        ,all_rat{i},...
       60,'filled','MarkerEdgeColor',fc(i,:),...
       'MarkerFaceColor',fc(i,:));
    hold on
end
xticklabels([]);

l2 = legend(pl,all_met_names,'position',[0.31 0.13 0.1 0.1],'fontsize',18);

pause(1); % DO NOT DELETE
for k = 1:length(l2.EntryContainer.NodeChildren)
    l2.EntryContainer.NodeChildren(k).Icon.Transform.Children.Children.Size = 12;
end
legend boxoff
ylabel({'Relative #', 'of electrodes','in 95% CI'})
%ylim([1,20])
set(gca,'fontsize',20)



%% Next plot is summary stats for global metrics)
axes(ha(11))
pl = zeros(size(all_95_ci_width,2),1);
for i = 1:size(all_95_ci_width,2)
    pl(i) = scatter(i*ones(size(all_95_ci_width,1),1)+0.05*randn(size(all_95_ci_width,1),1)...
        ,all_95_ci_width(:,i),...
       60,'filled','MarkerEdgeColor',fc(i+5,:),...
       'MarkerFaceColor',fc(i+5,:));
    hold on
    
end
xticklabels([]);
l3 = legend(pl,global_names_all,'location','northeast','fontsize',18);
legend boxoff
pause(1) % do not delete
for k = 1:length(l3.EntryContainer.NodeChildren)
    l3.EntryContainer.NodeChildren(k).Icon.Transform.Children.Children.Size = 12;
end
ylabel({'95% CI width', 'for metric'})
set(gca,'fontsize',20)


annotation('textbox',[0 0.88 0.1 0.1],'String',...
    'A','FontSize',35,'linestyle','none');
annotation('textbox',[0 0.63 0.1 0.1],'String',...
    'B','FontSize',35,'linestyle','none');
annotation('textbox',[0 0.38 0.1 0.1],'String',...
    'C','FontSize',35,'linestyle','none');


annotation('textbox',[0 0.13 0.1 0.1],'String',...
'D','FontSize',35,'linestyle','none');
annotation('textbox',[0.59 0.13 0.1 0.1],'String',...
'E','FontSize',35,'linestyle','none');


%print(gcf,[outFolder,'all_',freq,contig_text,sec_text],'-depsc');
%print(gcf,[outFolder,'all_',freq,contig_text,sec_text],'-dpng');

end