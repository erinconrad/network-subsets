function y = get_perc_elecs(X,n)

% This function finds the electrodes y that account for the nth percentile
% of all electrodes in the set X

n = round(n/100*length(X));

% get the unique elements of X and the times they first occur
X = sort(X);
[A,ia] = unique(X);

% Get the number of times each appears
num_each = [diff(ia);length(X)-ia(end)+1];

% Sort by how often they appear
[num_each,I] = sort(num_each,'descend');
most_common = A(I);

% Get the cumulative sum;
cs = cumsum(num_each);

% Get the index when the cumulative sum is >=n
greater = find(cs >= n);
first_greater = greater(1);

y = most_common(1:first_greater);

end