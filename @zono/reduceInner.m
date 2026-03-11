function obj = reduceInner(obj,final_order)
% REDUCEINNER Inner approximation order reduction for a zonotope

    G = obj.G;
    c = obj.c;

    n_g = size(G,2);
    n_r = final_order;

    %% Sort generators by magnitude
    norms = vecnorm(G);
    [~,indx_sort] = sort(norms,'descend');
    G = G(:,indx_sort);

    %% Split generators
    G1 = G(:,1:n_r);
    G2 = G(:,n_r+1:end);

    %% Compute dot products
    alpha = G1' * G2;
    alpha_abs = abs(alpha);

    %% Normalize dot products
    alpha_norm = alpha .* (1 ./ max(alpha_abs,[],1));

    %% Construct transformation matrix
    T2 = zeros(n_r,n_g-n_r);
    T2(alpha_norm == 1) = 1;
    T2(alpha_norm == -1) = -1;

    T = [eye(n_r); T2'];

    %% Transform generators
    G_transformed = G * T;

    %% Create reduced zonotope
    obj = zono(G_transformed,c);

end