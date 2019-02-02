function SMC = simple_matching_coefficient(A,B)

%{
If A and B are binary arrays, then the SMC is the number of matching
attributes divided by the number of attributes.

%}

if length(A) ~= length(B)
    error('Need A and B to be of equal length\n');
end

if strcmp(class(A),'logical') == 0 || strcmp(class(B),'logical') == 0
    error('A and B must be binary\n');
end

% Subtract A from B
diffAB = A-B;

% Get number of elements that are 0 (meaning that A and B matched at those
% positions)
n_equal = sum(diffAB == 0);

% Get SMC
SMC = n_equal/length(A);

end