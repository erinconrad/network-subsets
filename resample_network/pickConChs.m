function out_chs = pickConChs(locs,n,plotStuff,jitter,i)

%{
The goal of this function is to select 
contiguous electrodes. It takes the indicated channel i and returns its n
nearest neigbors.

Jitter is 0 and so there is no randomness.

%}

%% Parameters

chs = 1:size(locs,1);

if n > size(locs,1)
    error('Requesting too many electrodes\n');
end


% We will systematically loop through all channels
ch = i;

%{
    % Pick a random channel
    ch = randi(length(chs));
%}

    
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