function nodal_CI(stats,pt)

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

end