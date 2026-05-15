
function Z_reduced = reduceOuter(Zi,cons_to_remove, varargin)

    % Parse optional parameters
    p = inputParser;
    addParameter(p, 'tol', 1e-15, @isnumeric);
    addParameter(p, 'tolNonZero', 1e-6, @isnumeric);
    %tolr 1e-5 gave wrong results for several leaves. example leaf 13 
    addParameter(p, 'tolR', 1e-8, @isnumeric);
    addParameter(p, 'tolFeas', 1e-4, @isnumeric);
    addParameter(p, 'tolE', 1e-15, @isnumeric);
    addParameter(p, 'plotResult', false, @islogical);
    addParameter(p, 'plotComparison', false, @islogical);
    parse(p, varargin{:});
    
    % Initialize variables
    G = Zi.G;
    c = Zi.c;
    A = Zi.A;
    b = Zi.b;
    tol = p.Results.tol;
    tolNonZero = p.Results.tolNonZero;
    tolR = p.Results.tolR;
    tolFeas = p.Results.tolFeas;
    tolE = p.Results.tolE;
    
    % Display initial counts
    fprintf('%d generators and %d constraints (initial set).\n', size(G,2), size(A,1));

    % [A, b, G, ~] = refined_rref_with_pivoting(A, b, G);
    
    % Main reduction loop
    % while true
    if size(Zi.A,1)<cons_to_remove
        disp('Cannot remove more constaints than it exists')
        Z_reduced=Zi;
        return
    end

    for cons =1:cons_to_remove
        
        % Remove zero columns/rows
        zeroCols = find(all(abs([G; A]) < tol, 1));
        G(:, zeroCols) = [];
        A(:, zeroCols) = [];
        
        zeroRows = find(all(abs([A, b]) < tol, 2));
        A(zeroRows, :) = [];
        b(zeroRows) = [];
        
        % Calculate bounds using LP
        E = bound_e_LP_calculator(G, A, b);
        
        % Remove factors with equal bounds
        equal_indices = find(E(:,1) == E(:,2));
        if ~isempty(equal_indices)
            c = c + G(:, equal_indices) * E(equal_indices, 1);
            G(:, equal_indices) = [];
            b = b - A(:, equal_indices) * E(equal_indices, 1);
            A(:, equal_indices) = [];
            % continue; % Restart loop after removing equal bounds
        end
        
        % Calculate R bounds using LP
        R = bound_r_LP_calculator(G, A, b);
        
        % % Check for redundant factors
        % eff_r = find((R(:,1) >= -1 + tolR) & (R(:,2) <= 1 - tolR));
        % if isempty(eff_r)
        %     break; % No more redundancies found
        % end
        % 
        % Remove the last redundant factor found
        [~,order] = sort(max(abs(A),[],2),'descend');
        A = A(order,:);
        b = b(order);
        candidates = find(A(1,:) ~= 0);


        remove_position = chooseFactorToRemove(G, R(:,1), R(:,2), candidates);

        % remove_position = eff_r(end);
        
        
        % if isempty(removable_constraint)
        %     continue;
        % end
        % removable_constraint = removable_constraint(1);
        
        removable_constraint=1;

        % Remove the generator and constraint
        [G, A, b, c] = removeGenerator(G, A, b, c, remove_position, removable_constraint);
    % end
    end

    % original_dimension=size(G,1);
    % original_constrant_number=size(A,1);
    % 
    % obj=zono([G;A],[c;-b]);
    % 
    % obj= outerApproximation(obj);
    % 
    % [G,c,A,b] = splitLiftedZonotope(obj.G, obj.c, original_dimension, original_constrant_number);

    % Gtilde=obj.G;
    % ctilde=obj.c;
    % 
    % G=Gtilde(1:original_dimension,:);
    % A=Gtilde(original_dimension+1:end,:);
    % 
    % c=ctilde(1:original_dimension);
    % b=ctilde(original_dimension+1:end);

    % Create output zonotope
    Z_reduced = conZono(G, c, A, b);
        
end


function [G_new, A_new, b_new, c_new] = removeGenerator(G, A, b, c, remove_position, removable_constraint)
%REMOVEGENERATOR Remove redundant generator and corresponding constraint
%
%   Inputs:
%       G - Generator matrix
%       A - Constraint matrix
%       b - Constraint vector
%       c - Center vector
%       remove_position - Index of generator to remove
%       removable_constraint - Index of constraint to remove
%
%   Outputs:
%       G_new - Updated generator matrix
%       A_new - Updated constraint matrix
%       b_new - Updated constraint vector
%       c_new - Updated center vector

    [nc, ng] = size(A); % Number of generators and constraints
    
    % Step 1: Define the E matrix
    E_lam = zeros(ng, nc);
    E_lam(remove_position, removable_constraint) = 1; % Place a 1 at the specified position
    
    % Step 2: Compute Lambda matrices
    a_1j = A(removable_constraint, remove_position); % The pivot element
    Lambda_G = G * E_lam / a_1j; % Lambda_G = G * E * a_1j^(-1)
    Lambda_A = A * E_lam / a_1j; % Lambda_A = A * E * a_1j^(-1)
    
    % Step 3: Update the generator matrix G
    G_new = G - Lambda_G * A;
    
    % Step 4: Update the center vector c
    c_new = c + Lambda_G * b;
    
    % Step 5: Update the constraint matrix A
    A_new = A - Lambda_A * A;
    
    % Step 6: Update the constraint vector b
    b_new = b - Lambda_A * b;
    
    % Step 7: Remove the selected generator and constraint
    G_new(:, remove_position) = []; % Remove generator column
    A_new(removable_constraint, :) = []; % Remove constraint row
    A_new(:, remove_position) = []; % Remove generator column from A
    b_new(removable_constraint) = []; % Remove constraint element
end



function j_best = chooseFactorToRemove(G, rhoL, rhoU, candidates)
% chooseFactorToRemove
%
% Selects the best generator/factor to remove during constrained zonotope
% constraint reduction using the bound-based heuristic from Scott (2016).
%
% INPUTS
% G           : generator matrix (n x ng)
% rhoL        : lower bounds on factors (ng x 1)
% rhoU        : upper bounds on factors (ng x 1)
% candidates  : indices of factors that can be removed
%
% OUTPUT
% j_best      : index of the best factor to eliminate

ng = size(G,2);

bestScore = inf;
j_best = candidates(1);

for k = 1:length(candidates)

    j = candidates(k);

    % amount by which constraint forces factor outside [-1,1]
    rj = max(0 , max(abs(rhoL(j)),abs(rhoU(j))) - 1);

    % if zero error possible choose immediately
    if rj == 0
        j_best = j;
        return
    end

    % generator influence
    gj = norm(G(:,j),2);

    % error estimate (Hausdorff approximation)
    score = gj * rj;

    if score < bestScore
        bestScore = score;
        j_best = j;
    end

end
end


% function obj= outerApproximation(obj)
%     G = obj.G;
%     c=obj.c;  
%     X = [G, -G]';  
%     Co = X' * X;   
%     % Step 4: Perform Singular Value Decomposition (SVD)
%     [U, ~, ~] = svd(Co);  % U contains the principal components (eigenvectors)    
%     Z_transformed=U'*obj;    
%     obj=U*boundingBox(Z_transformed);
% end
% 
% 
% function [Gnew,cnew,Anew,bnew] = splitLiftedZonotope(Gplus_red, cplus_red, n, nc)
% % Gplus_red : reduced lifted generator matrix, size (n+nc) x q
% % cplus_red : reduced lifted center, size (n+nc) x 1
% % n         : original state dimension
% % nc        : number of constraint rows kept in lifted form
% 
% cnew = cplus_red(1:n);
% bnew = -cplus_red(n+1:n+nc);
% 
% Gnew = Gplus_red(1:n,:);
% Anew = Gplus_red(n+1:n+nc,:);
% 
% end