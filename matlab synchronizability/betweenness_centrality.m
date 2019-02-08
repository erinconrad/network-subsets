function bc = betweenness_centrality(A,method)

%{
Need to make sure I am interpreting edge weight correctly

https://www.sciencedirect.com/science/article/pii/S1388245718302451?via%3Dihub
Betweenness centrality of intracranial electroencephalography networks and surgical epilepsy outcome

%}


%% Method 1: Weighted graph
G = graph(A);
bc1 = centrality(G,'betweenness','Cost',1./G.Edges.Weight);

%% Method 2: Threshold
thresh_perc = 0.8;
con_thresh = quantile(A(:),thresh_perc);
A(A<con_thresh) = 0;
A(A>con_thresh) = 1;
G = graph(A);
bc2 = centrality(G,'betweenness');

if method == 1
    bc = bc1;
elseif method == 2
    bc = bc2;
end

%{
scatter(bc1,bc2)
corr(bc1,bc2)
%}

end