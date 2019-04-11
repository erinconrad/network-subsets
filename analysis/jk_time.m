function jk_time(stats,pt)

%{

This creates the supplemental figure looking at the changes in jackknife
results for different times

%}

%% Parameters
contig_text = 'random';
freq = 'high_gamma';
all_sec = {'sec_neg10','sec_neg5','sec_0','sec_5','sec_10'};
sec_texts = {'EEC - 10s','EEC - 5s','EEC','EEC + 5s','EEC + 10s'};
which_texts = {'ns','min_cc_elecs'};
which_names = {'Node strength','Regional control centrality'};
whichPt = 1;
which_global = {'sync','eff','trans'};
global_names_all = {'Synchronizability','Global efficiency','Transitivity'};

%% Locations

[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
outFolder = [resultsFolder,'cc_comparison/'];
locs = pt(whichPt).new_elecs.locs;

figure
set(gcf,'Position',[50 0 1400 900]);
[ha,pos] = tight_subplot(3, 5, [0 0.01], [0.05 0.04],[0.045 0.02]);
set(ha(11),'position',[pos{11}(1), pos{11}(2), ...
    (pos{5}(1)+pos{5}(3)-pos{1}(1))/2-0.03, pos{11}(4)-0.02]);
set(ha(12),'position',[pos{11}(1)+(pos{5}(1)+pos{5}(3)-pos{1}(1))/2+0.03,...
    pos{11}(2),(pos{5}(1)+pos{5}(3)-pos{1}(1))/2-0.03, pos{11}(4)-0.02]);
delete(ha(13))
delete(ha(14))
delete(ha(15))

count = 0;

%% Plot the first two rows
for text = 1:length(which_texts)
    for i = 1:length(all_sec)
        count = count + 1;
        base = stats(whichPt).(freq).(contig_text).(all_sec{i});
        axes(ha(count))
        new_count = 0;
        pl = zeros(5,1);
        c = colormap(parula(5));
        ex_single = [];
        ex_regional = [];
        if strcmp(which_texts{text},'min_cc_elecs') == 0
            true1 = base.(which_texts{text}).true';
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
            
            if strcmp(which_texts{text},'min_cc_elecs') == 0
                elecs = base.(which_texts{text}).(sprintf('single_%d',j))';
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
        
        
        scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k','linewidth',2);
        if count == 1 || count == 6
            zlabel(sprintf('%s',which_names{text}));
        elseif count == 10
           l1= legend(pl,{'True','70%','80%','90%','95%'},'Position',[0.94 0.58 0.03 0.10],...
            'box','on');
            pause(1)
            for k = 1:length(l1.EntryContainer.NodeChildren)
                l1.EntryContainer.NodeChildren(k).Icon.Transform.Children.Children.Size = 14;
            end
        end
        
        if ismember(count,[1:5]) == 1
            title(sprintf('%s',sec_texts{count}));
        end
        
        
                
        view(150,32.4);
        set(gca,'fontsize',20)
        xticklabels([])
        yticklabels([])
        zticklabels([])

    end
end

%% Plot the last row first column
for gl_idx = 1:2
    count = count + 1;
    axes(ha(count));

    % Get global for each
    true = zeros(length(all_sec),1);
    all_perm = zeros(length(all_sec),100);
    for i = 1:length(all_sec)

        base = stats(whichPt).(freq).(contig_text).(all_sec{i});
        true(i) = base.(which_global{gl_idx}).true;
        all_perm(i,:) = base.(which_global{gl_idx}).all;

    end

    violin(all_perm',true,'facecolor',[0 0.4470 0.7410])
    xticks(1:length(sec_texts))
    xticklabels(sec_texts)
    ylabel(sprintf('%s',global_names_all{gl_idx}));
    set(gca,'fontsize',20)
end

%% Add annotations
annotation('textbox',[0 0.87 0.1 0.1],'String',...
    'A','FontSize',35,'linestyle','none');
annotation('textbox',[0 0.55 0.1 0.1],'String',...
    'B','FontSize',35,'linestyle','none');
annotation('textbox',[0 0.25 0.1 0.1],'String',...
    'C','FontSize',35,'linestyle','none');
annotation('textbox',[0.5 0.25 0.1 0.1],'String',...
    'D','FontSize',35,'linestyle','none');

print(gcf,[outFolder,'time_sens'],'-depsc');
print(gcf,[outFolder,'time_sens'],'-dpng');

end