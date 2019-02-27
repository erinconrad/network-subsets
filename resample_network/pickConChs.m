function out_chs = pickConChs(locs,n,plotStuff,jitter)

%{
The goal of this function is to select a random group of n more or less
contiguous electrodes. It picks a random channel, and then sorts all
channels based on their distance from this channel, and then it selects the
closest n of them. It also adds some random jitter to the distances so that
it gets random groups of electrodes.

%}

%% Parameters
% jitter: 0-10 nearly contiguous; 100+ very discontiguous

chs = 1:size(locs,1);

if n > size(locs,1)
    error('Requesting too many electrodes\n');
end


% Pick a random channel
ch = randi(length(chs));

% Get distances to other chs
dist = vecnorm(locs-repmat(locs(ch,:),length(chs),1),2,2);

% Random jitter to lie about distances
x = randi([-1,1],length(dist),1)*jitter;
dist = dist + x;

% Re-sort channels by distances to ch
[~,I] = sort(dist);
chs_sorted = chs(I);
    
% Take the first n channels
out_chs = chs_sorted(1:n);

if plotStuff == 1
    figure
    scatter3(locs(:,1),locs(:,2),locs(:,3),100,'k');
    hold on
    scatter3(locs(out_chs,1),locs(out_chs,2),locs(out_chs,3),100,[1 0 0],'filled');
    pause
    close(gcf)
    
end


end