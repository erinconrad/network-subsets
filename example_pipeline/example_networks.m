function example_networks

%% Need to edit

% Edit this to match the file path
data_folder = '../../data/example_data/';

% Edit this with the path to the Brain Connectivity Toolbox (necessary to run this code)
BCT_folder = '../../BCT/';

%% Other parameters
% Change if desired
which_times = 0; %0 is EEC; can also choose from [-10, -5, 0, 5, 10] (how many seconds after EEC)
freq = 'high_gamma'; % Can also do 'beta'

% Do not change
adj_file = 'HUP083_ex.mat'; % pt 11
info_file = 'info.mat';
whichPt = 11;
sz_num = 7;

%% Add script paths
addpath(genpath(BCT_folder)); % add the BCT folder
addpath(genpath('./..')); % add the main script folder

%% Load files
pt = load([data_folder,info_file]);
pt = pt.pt;

adj = load([data_folder,adj_file]);
adj = adj.adj;

%% Perform network resampling to get reliability of nodal metrics
% We can't obtain global metric reliability with just one patient

% Prep data
example.adj = adj;
example.pt = pt;
example.whichPt = whichPt;
example.sz_num = sz_num;
example.which_times = which_times;
example.freq = freq;

% Call main network stats function to do the main resampling
out = network_stats([],0,[],[],example); 

% Analysis
compare_metrics(pt,out(whichPt),1);

%% Get confidence intervals of global and nodal metrics
fprintf(['Showing confidence intervals for single patient of:\n'...
    'A: Location of highest node strength electrodes\n'...
    'B: Location of lowest regional control centrality electrodes\n'...
    'C: 95%% confidence interval for synchronizability\n'...
    'D: Number of electrodes forming 95%% CI set of nodal metrics\n'...
    'E: 95%% confidence interval for all global metrics\n']);
all_CI(out(11),pt,1)

end