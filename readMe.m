%% Network subsets analysis
%{ 
Written by Erin Conrad and John Bernabei, 2019, University of
Pennsylvania.

This codebase contains the scripts for performing the analysis associated
with the paper "The sensitivity of network statistics to incomplete 
electrode sampling on intracranial EEG". 

The information below explains the components of the codebase.

To run an example script, navigate to the folder example_pipeline, open the
script example_networks, and read the instructions. This is an example
script demonstrating how the codebase is run using example data from a
single patient and single seizure.
 
 %}

%% Overview
%{

This analysis determines the sensitivity/robustness of network theory
results to incomplete spatial sampling generated by sparse electrode
placement.

To run this, you will need adjacency matrices that represent the functional
connectivity of electrodes (coherence for frequency bands of interest).
These take the form of undirected, weighted, symmetric adjacency matrices.
Once you have these, do the following steps:

1) Make pt struct (or take our pre-made pt structure "info.mat", provided
as part of the example_script).
2) Convert python adjacency matrices to .mat files (see correct .mat file
format in the example "HUP083_ex.mat", provided as part of the
example_script)
3) Get electrode locations aligned with adjacency matrix format
4) Do resampling and get network changes
5) Do high-level analyses

%}

%% 1: Make pt struct
%{
Go to struct_setup/ and run:

- getPtStruct: this is a one time use function that takes electrode data
and a clinical json file (need to have these available) and creates a
patient structure. The electrode identities will NOT be appropriately
aligned with the adjacency matrix electrode identities.

%}

%% 2: Convert python adjacency matrices to .mat files
%{
The adjacency matrices were originally in .npz format. They will need to be
.mat files to use in this pipeline. To convert them, run the following
script in convertPython/:

- unzip_npz_files: takes a patient structure and a list of which patients to run,
and generates .mat files containing the electrode info and adjacency
matrices from the .npz file
%}

%% 3: Get electrode locations aligned with adjacency matrix format
%{
Once you have a patient structure and adjacency matrices as .mat files, you
need to redefine the electrode identities in the patient structure so that
they line up with the adjacency matrices. Go to struct_setup and run:

- redefineElecs (requires the patient structure and the adjacency matrices)

%}

%% 4: Do resampling and get network changes
%{
This is the bulk of the run time. This is what actually goes through,
calculates network metrics, removes electrodes, recalculates the metrics,
and compares them with the original metrics. Navigate to networkStats/ and
run:

- network_stats: specify which patients to do it for (1:33 is all), whether
you are making the main structure (do_soz_analysis = 0), making the
structure to correlate metric agreement with distance from resection zone
(do_soz_analysis = 1), making the structure to compare metric agreement for 
SOZ-sparing vs SOZ-targeted resampling (do_soz_analysis = 2), what seizure 
you are doing (1, 2, or 3, representing the first, second, or last seizure,
respectively), what network density you want, what time window you want,
and finally, whether you are running the code as part of the example
script.

This has many dependencies, the most important of which is:

- resampleNetwork (in resample_network/): this takes the adjacency matrix
and info about number of electrodes to remove and number of iterations to
do this over, resamples the network according to these rules, and
calculates new network metrics

- There are three potential outputs to network_stats, depending on whether
do_soz_analysis is 0, 1, or 2:
- stats.mat: this is the output if do_soz_analysis is 0. This contains info
on metric reliablity as well as the result of the jackknife analysis
(generating confidence intervals).
- soz.mat: this is the output if do_soz_analysis is 1. This contains metric
agreement as a function of distance of the removed electrodes from the
resection zone.
- soz_overlap.mat: this is the output if do_soz_analysis is 2. This
compares metric agreement between seizure onset zone-targeting and seizure
onset zone-sparing subsampling.

%}

%% 5: Do high-level analyses
%{
There are three major analyses, all of which can be run by going to
analysis/:

- compare_metrics: this takes the stats.mat structure and compares (with
stats and with plots) the reliability of the different metrics
- compare_soz_resec: this takes the soz.mat structure and assesses (with
stats and with plots) the correlation between the metric agreement and the
distance from the seizure onset zone
- soz_overlap_analysis: this takes the soz_overlap.mat structure and
assesses whether resampling that targets the SOZ has a larger effect on the
network than resampling that spares the SOZ
- all_CI: this takes the stats.mat structure and generates plots exploring
the jackknife approach 
- stats_CI: this takes the stats.mat structure and generates stats
summarizing the results of the jackknife approach

% More minor functions in the analysis/ folder:
- basic_info_networks: gives basic info about the patients
- jk_time: gets info about how the jackknife results change across time.

%}





