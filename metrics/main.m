clear
rng default

%% Parameters
freq = [100 150];
coherence_window_length = 100;
coherence_sample_overlap = 80;
Fs = 1000;

%% Load in data

% This is all fake data
t = 0:1/Fs:1-1/Fs;
x = cos(2*pi*100*t)+sin(2*pi*200*t)+0.5*randn(size(t));
y = 0.5*cos(2*pi*100*t-pi/4)+0.35*sin(2*pi*200*t-pi/2)+ ...
    0.5*randn(size(t));
z = 0.8*cos(2*pi*100*t-pi/4)+0.65*sin(2*pi*200*t-pi/2)+ ...
    0.6*randn(size(t));
data = [x;y;z]';


%% Clean it up



%% Common average reference

% Subtract signal averaged over all electrodes
data = data - repmat(mean(data,2),1,size(data,2));



%% Calculate adjacency matric A using spectral coherence


A = spectral_coherence(data, ...
    freq,coherence_window_length,...
    coherence_sample_overlap,...
    Fs);



%% Calculate synchronizability

sync = synchronizability(A);
