function [A,E] = makeFakeA(N,p)

A = zeros(N,N);
E = [];
for i = 1:size(A,1)
    for j = 1:i-1
        makeOne = randi(100);
        if makeOne < p*100
            A(j,i) = 1;
            A(i,j) = 1;
            E = [E;j i];
        end
    end
end


end