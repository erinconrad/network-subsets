function make_heatmap(which_tbl)

%% Parameters
if which_tbl == 4 || which_tbl == 5 || which_tbl == 6
    color_range = [-4 4];
elseif which_tbl == 3
    color_range = [0.6 1];
elseif which_tbl == 2
    color_range = [0.7 1];
else
    color_range = nan;
end

%% Locations
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
tblFolder = [resultsFolder,'supplemental_tables/'];
toolFolder = [mainFolder,'tools/'];
addpath(genpath(toolFolder));

%% Load table
file = [tblFolder,'table_',sprintf('%d',which_tbl),'.xlsx'];
T = readtable(file,'ReadRowNames',true);

%% Get correct column names
column_names = cell(length(T.Properties.VariableDescriptions),1);
for i = 1:length(T.Properties.VariableDescriptions)
    temp_name = T.Properties.VariableDescriptions{i};
    new_name = erase(temp_name,'Original column heading: ');
    new_name = erase(new_name,'''');
    if isempty(new_name) == 1
        column_names{i} = T.Properties.VariableNames{i};
    else
        column_names{i} = new_name;
    end
end

%% Get matrix of just numbers and cell with asterixes

% Convert table to cell
C = table2cell(T);

% Initialize new cells with numbers and sig values
c_num = zeros(size(C));
c_sig = C;

% Loop through elements of cell
for i = 1:size(C,1)
    for j = 1:size(C,2)
        
        element = C{i,j};
        
        % Get class of element
        cl = class(element);
        
        
        if strcmp(cl,'double') == 1
        
            % It's not significant
            c_sig{i,j} = sprintf('%1.2f',element);
            c_num(i,j) = element;
                
        elseif strcmp(cl,'char') == 1
            
            % Get number asterixes
            num_stars = count(element,'*');
            if num_stars == 0
                c_sig{i,j} = sprintf('%1.2f',str2double(element));
                no_star = element;
            else
                
                no_star = erase(element,repmat('*',1,num_stars));
                c_sig{i,j} = sprintf('%1.2f%s',str2double(no_star),repmat('*',1,num_stars));
            end
            
            c_num(i,j) = str2double(no_star);
                   
        else
            
            error('what\n');
            
        end
        
    end
end

%% Make heatmap
if isnan(color_range) == 1
    h = heatmap_custom(c_num,column_names,T.Properties.RowNames,c_sig,...
    'fontsize',15,'TickFontSize',12,'ShowAllTicks',1);
else
    h = heatmap_custom(c_num,column_names,T.Properties.RowNames,c_sig,...
    'MinColorValue',color_range(1),'MaxColorValue',color_range(2),...
    'fontsize',15,'TickFontSize',12,'ShowAllTicks',1);
end

colorbar

if which_tbl == 6
    set(gcf,'position',[10 313 1390 485]);
elseif which_tbl == 4
    set(gcf,'position',[10 313 1350 485]);    
elseif which_tbl == 2
    set(gcf,'position',[10 313 1300 485]);   
elseif which_tbl == 5
    set(gcf,'position',[10 313 600 485]);  
elseif which_tbl == 3
    set(gcf,'position',[10 313 600 485]);  
end

%set(gca,'LooseInset',get(gca,'TightInset'));

% Save figure
print(gcf,[tblFolder,sprintf('heatmap_%d',which_tbl)],'-depsc');

end
