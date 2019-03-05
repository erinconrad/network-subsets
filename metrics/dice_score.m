function ds = dice_score(X,Y)

% https://en.wikipedia.org/wiki/S%C3%B8rensen%E2%80%93Dice_coefficient
% The dice score is twice the number of elements common to both sets 
% divided by the sum of the number of elements in each set

ds = 2*(sum(X==1 & Y ==1))/(sum(X) + sum(Y));


end
