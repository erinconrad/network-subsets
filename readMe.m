%% Network subsets analysis

%{
These functions take an adjacency matrix A, and randomly remove nodes (so
associated rows and columns from the matrix), and determines how the
control centralities of each node changes in the resampled network.

To run, navigate to networkStats and run network_stats, inputting the
patient number you want to run. This will require adjacency matrices and
electrode locations to be stored in their own structures in the data folder

%}