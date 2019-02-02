function c_c = control_centrality(A)

n_ch = size(A,1);

% Calculate the synchronizability
sync = synchronizability(A);

% To calculate control centrality of each electrode, lesion it out
% and re-calculate the sync
c_c = zeros(n_ch,1);
for ich = 1:n_ch

    % Remove the channel
    A_temp = A;
    A_temp(ich,:) = [];
    A_temp(:,ich) = [];

    % Recalculate the synchronizability
    sync_temp = synchronizability(A_temp);

    % Calculate control centrality
    c_c(ich) = (sync_temp - sync)/sync;

end

end