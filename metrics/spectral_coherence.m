function A = spectral_coherence(data,freq,window_length,sample_overlap,Fs)

n_chan = size(data,2);

% Initialize adjacency matrix
A = zeros(n_chan, n_chan);

% Compute all coherences
for i = 1:n_chan
   for j = 1:n_chan
      [out,F] =  mscohere(data(:,i),data(:,j),...
          hamming(window_length),...
          sample_overlap,...
          window_length,...
          Fs);
      
      % Find the indices of the frequency band closest to the desired frequency band
      
       closest_f = intersect(find(F>=freq(1)),find(F<=freq(2)));
       
       % Store coherence in A
       A(i,j) = mean(out(closest_f));
       
   end
    
end



end