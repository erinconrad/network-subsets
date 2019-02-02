function s_c = structConsistency(A)

%{
This function calculates the structural consistency of an undirected
binary network with adjacency matrix A. The algorithm is taken from:

Lü, Linyuan, et al. "Toward link predictability of complex networks." 
Proceedings of the National Academy of Sciences 112.8 (2015): 2325-2330.

%}

%% Parameters

% How many permuations to do
nboot = 1e2;

% What percentage of links to remove to make delta
pH = 0.1;

% Plot every new matrix? Just for QA purposes
doPlot = 0;

%% Generate fake adjacency matrix if needed
if isempty(A) == 1
    fprintf('Warning, no adjacency matrix given, will use fake data.\n');
    
    %{
    A = zeros(9,9);
    E = [1 2;1 4;2 3;2 5;4 5;3 6;5 6;4 7;...
        5 8;6 9;7 8;8 9];
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            for k = 1:size(E,1)
                if i == E(k,1) && j == E(k,2)
                    A(i,j) = 1;
                    A(j,i) = 1;
                end
            end
        end
    end
    %}
    
    %% Make fake A and E
    
    % Size of matrix
    N = 80;
    
    % Probability of link connection
    p = 0.2;
   
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

% How many links to take as delta E
n_links = ceil(size(E,1)*pH);

s_c_boot = zeros(nboot,1);
for ib = 1:nboot

    %% Select random perturbation set E_delta and generate A_r and A_d
    which_rows = randperm(size(E,1),n_links);
    E_d = E(which_rows,:);
    
    E_r = E(~ismember(E, E_d, 'rows'),:);
    A_r = A;
    A_d = zeros(size(A));

    for i = 1:size(A,1)
        for j = 1:size(A,2)
            for k = 1:size(E_d,1)
                if i == E_d(k,1) && j == E_d(k,2)
                    A_r(i,j) = 0;
                    A_r(j,i) = 0;
                    A_d(i,j) = 1;
                    A_d(j,i) = 1;
                end
            end
        end
    end

    %% Check that A = A_r + A_d
    if isequal(A_r+A_d,A) == 0
        error('What\n');
    end

    %% Calculate the eigenvalues and eigenvectors of A_r
    [x_k,l_k,~] = eig(A_r);
    l_k_nums = eig(A_r);

    %% Check that x_k and l_k make sense as eigenvalues and eigenvectors of A_r
    if sum(sum(abs(A_r*x_k-x_k*l_k))) > 1e-5
        error('what\n');
    end

    %% Check equation 1 in paper
    x_k_t = x_k.';
    tot_sum = 0;
    for k = 1:length(l_k_nums)
        tot_sum = tot_sum + l_k_nums(k)*x_k(:,k)*x_k_t(k,:);
    end

    if abs(A_r-tot_sum) > 1e-5
        error('What\n');
    end

    %% Calculate l_d
    l_d = zeros(size(l_k_nums));
    for k = 1:length(l_k_nums)
        %l_d = x_k.'*A_d*x_k/(x_k.'*x_k);
        l_d(k) = x_k_t(k,:)*A_d*x_k(:,k)/(x_k_t(k,:)*x_k(:,k));
    end

    %% Calculate perturbed matrix
    A_p = 0;
    for k = 1:length(l_k_nums)
        A_p = A_p + ...
            (l_k_nums(k) + l_d(k))*x_k(:,k)*x_k_t(k,:);
    end

    %% Make near-zero values zero
    for i = 1:size(A_p,1)
        for j = 1:size(A_p,2)
            if abs(A_p(i,j)) < 1e-5
                A_p(i,j) = 0;
            end
        end
    end


    %% Construct U, universal link array, with A_p scores
    U = [];
    A_p_scores = [];
    count = 0;
    for i = 1:size(A,1)
        for j = 1:i-1
            count = count + 1;
            U = [U;j i];
            A_p_scores = [A_p_scores;count A_p(i,j)];
        end
    end

    if size(U,1) ~= (size(A,1)^2-size(A))/2
        error('What\n');
    end

    %% Remove E_r
    E_no = U(~ismember(U, E_r, 'rows'),:);
    A_p_scores_no = A_p_scores(~ismember(U, E_r, 'rows'),:);

    %% Rank remaining links by A_p scores
    [~,I] = sort(A_p_scores_no(:,2),'descend');
    E_no_sorted = E_no(I,:);

    %% Get first L links, where L is the number of links in E_d
    E_L = E_no_sorted(1:size(E_d,1),:);

    %% Show the locations of these links
    if doPlot == 1
    figure
    set(gcf,'Position',[136 567 1305 231]);
    subplot(1,4,1)
    imagesc((A))
    colormap(flipud(gray))
    title('Original matrix');

    subplot(1,4,2)
    imagesc((A_d))
    colormap(flipud(gray))
    title('Removed links');

    subplot(1,4,3)
    imagesc(A_p)
    colormap(flipud(gray))
    title('Perturbed matrix')

    A_show = zeros(size(A));
    for i = 1:size(A,1)
        for j = 1:size(A,2)
            if ismember([i,j],E_L,'rows') == 1
                A_show(i,j) = 1;
                A_show(j,i) = 1;
            end
        end
    end

    subplot(1,4,4)
    imagesc(A_show)
    colormap(flipud(gray))
    title('Highest ranked links');
    
    pause
    close(gcf)
    end


    %% Define structural consistency
    s_c_boot(ib) = sum(ismember(E_d,E_L,'rows'))/size(E_d,1);
end
    
s_c = mean(s_c_boot);

if 1 == 0
    figure
    scatter(1:nboot,sort(s_c_boot));
end

fprintf('Structural consistency %1.1f.\n',s_c);

end