function example_networks

%{
This function executes some of the primary analyses in the paper "The 
sensitivity of network statistics to incomplete electrode sampling on 
intracranial EEG". The goal is to provide an example of the structure of
the input data and how the code utilizes this data. 

To run this code, you will need the following things:
- the "Network Subsets" codebase, available at https://github.com/erinconrad/network-subsets/
- The Brain Connectivity Toolbox (Complex network measures of brain 
connectivity: Uses and interpretations. Rubinov M, Sporns O (2010) 
NeuroImage 52:1059-69.), available at https://sites.google.com/site/bctnet/
- the patient information "info.mat" file, which contains de-identified
information about patients. This is available for download along with the
primary codebase on Github in the example_pipeline folder.
- the example adjacency matrix file "HUP083_ex.mat", which contains the
adjacency matrices for a single seizure from a single patient (HUP083),
available along with the code base in the example_pipeline folder.

Once you have downloaded all of these components, edit the "Paths to edit" 
section below to point to where you have saved the "info.mat" file, the
"HUP083_ex.mat" file, and the Brain Connectivity Toolbox. Then navigate to
the folder containing this function and run
>> example_networks

Because this code only performs analyses for a single patient, it will not
perform across-patient level analyses. For instance, it will not calculate
global reliability metrics, which require calculating variance across 
patients.
%}

%% Paths to edit
% Edit this to match the file path to "info.mat" and "HUP083_ex.mat"
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
fprintf('Preparing data for analysis...\n\n');
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
fprintf('Performing main resampling method (this is the longest step)...\n\n');
out = network_stats([],0,[],[],[],example); 
fprintf('Got main network statistics. Press any key to continue to analysis...\n\n');

% Analysis to compare nodal metrics for this single patient
fprintf('Now comparing nodal metric reliability for this single patient...\n\n');
compare_metrics(pt,out(whichPt),1);
fprintf('Press any key to continue to jackknife resampling analysis...\n\n');


%% Get confidence intervals of global and nodal metrics
fprintf('Now showing jackknife resampling analysis for individual patient confidence intervals...\n\n')
all_CI(out(whichPt),pt,1)
fprintf(['Showing confidence intervals for single patient of:\n'...
    'A: Location of highest node strength electrodes\n'...
    'B: Location of lowest regional control centrality electrodes\n'...
    'C: 95%% confidence interval for synchronizability\n'...
    'D: Number of electrodes forming 95%% CI set of nodal metrics\n'...
    'E: 95%% confidence interval for all global metrics\n\n'...
    'Press any key to continue to the seizure onset zone analyses...\n\n']);

%% Perform network resampling to do seizure onset zone analyses
% Call main network stats function to do the main resampling
fprintf('Now performing network resampling for seizure onset zone analyses...\n\n');
out_soz_1 = network_stats([],1,[],[],[],example); % This does the distance-agreement association resampling
out_soz_2 = network_stats([],2,[],[],[],example); % This does the SOZ-sparing vs SOZ-targeted resamping
fprintf('Finished resampling. Press any key to analyze results...\n\n');

% Analysis to see if agreement between original and resampled metric
% correlates with distance from seizure onset zone
fprintf(['Now doing the analysis correlating the agreement between the\n'...
    'original and resampled metric with the distance of the removed electrodes\n'...
    'from the seizure onset zone...\n\n']);
compare_soz_resec(out_soz_1(whichPt),pt,1,0)
fprintf(['This figure shows the correlation between the resampled-original\n'...
    'metric agreement and the distance of the removed electrode contacts\n'...
    'from the seizure onset zone.\n\nPress any key to continue to the next analysis...\n']);
pause


% Analysis to see if agreement between original and resampled metric is
% better when we spare the seizure onset zone from removal
fprintf(['Now doing the analysis comparing seizure onset zone-sparing to\n'...
    'seizure onset zone-targeted resampling....\n\n']);
soz_overlap_analysis(out_soz_2(whichPt),pt,1,0,0)
fprintf(['This figure shows the agreement between the resampled and original'...
    '\nmetric for both seizure onset zone-targeted and seizure onset zone-sparing\n'...
    'resampling. For nodal metrics the agreement metric is the Spearman rank\n'...
    'correlation between the original and resampled metric. For global metrics it\n'...
    'is the negative of the absolute value of the relative difference between\n'...
    'the original and resampled metric.\n\nThis concludes the example pipeline.\n']);


end