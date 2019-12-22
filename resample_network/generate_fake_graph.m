function C = generate_fake_graph(A)

%{
This function generates a fake graph of the same size as A
%}

% Get the upper triangular matrix from A
U = triu(A,1);

% Get the number of elements in the upper triangle, excluding the diagonal
n_non_diag = (size(A,1)^2-length(A))/2; 

% Get the non-diagonal elemtns
C = logical(triu(ones(size(A)),1));
non_diag = A(C);

if length(non_diag) ~= n_non_diag
    error('What\n');
end

% Generate random permutation
p = randperm(n_non_diag);


% Generate a lower triangular matrix by permuting things
B = zeros(size(A));

count = 0;
for i = 1:size(U,1)
    for j = 1:i-1
        count = count + 1;
        B(i,j) = non_diag(p(count));
    end
end

if count ~= n_non_diag
    error('what\n');
end

[n,m] = size(B);
C = B + B';

end