function A = fake_A_certain_size(N)

A = zeros(N,N);

for i = 1:size(A,1)
    for j = 1:i-1
        a = abs(randn);
        A(i,j) = a;
        A(j,i) = a;
    end
end

end