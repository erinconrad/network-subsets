function test_transitivity(A)

% If I make a truly fake graph, transitivity doesn't change with graph size

% What if I make a more connected graph

%% Load Stuff
[electrodeFolder,jsonfile,scriptFolder,resultsFolder,...
pwfile,dataFolder,bctFolder,mainFolder] = resectFileLocs;
p1 = genpath(scriptFolder);
addpath(p1);

% Add brain connectivity toolbox
addpath([bctFolder]);

% Make a fake graph
%A = fake_A_certain_size(80);

% Get the transitivity
trans = transitivity_wu(A);

% Loop over 1,000 permutation
nb = 1000;
trans_fake = zeros(nb,1);
for i = 1:nb
    
    % Pick 20 random "electrodes"
    rm_idx = randsample(length(A),20);
    
    % remove them
    A_fake = A;
    A_fake(:,rm_idx) = [];
    A_fake(rm_idx,:) = [];
    
    % recalculate transitivity
    trans_fake(i) = transitivity_wu(A_fake);
    
end

% Plot distribution of transitivities
figure
plot(trans_fake,'o')
hold on
xl = get(gca,'xlim');
plot([xl(1) xl(2)],[trans trans])

mean(trans_fake)
trans

end