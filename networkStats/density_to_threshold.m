function threshold = density_to_threshold(density,A)

%{
This function takes a desired graph density, a graph A, and returns a
weighting threshold that will return that density
%}

% Graph size
N = size(A,1);

% Get upper triangular elements
ut = triu(A,1);

% Sort them; remove the zeros
ut_sort = sort(ut(:));
ut_sort = ut_sort(ut_sort > 0);

% Get the number of elements needed to be supra-threshold to achieve the desired density
n_needed = floor(length(ut_sort)-density*length(ut_sort));

% Make sure it's bounded by the number of elements
n_needed = min(n_needed, length(ut_sort));
n_needed = max(n_needed,1);


% Find the weight of the n_needed element of the array; this is our
% threshold
threshold = ut_sort(n_needed);

% Recalculate density (per BCT code)
A_new = A;
A_new(A_new <= threshold) = 0;
[dens_confirm,~,~] = density_und(A_new);
if abs(dens_confirm-density) > 0.001
    error('what\n');
end

% For testing, plot the new graph
if 0
    
    
    
    
    figure
    imagesc(A_new)
    colorbar
    title(sprintf('Density %1.2f, threshold %1.2f',dens_confirm,threshold))
    
end

end