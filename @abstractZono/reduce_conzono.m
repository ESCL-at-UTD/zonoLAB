function Z_reduced = reduce_conzono(Zi)
%REDUCECONZONOGUROBI Reduce redundant generators and constraints using Gurobi
% Reduction in 2D and 3D should most definitely work. Reduction in higher
% dimensions is still in question.
%
%   Input:
%       Zi - Constrained zonotope with fields: c, G, A, b
%
%   Output:
%       Z_reduced - Reduced constrained zonotope

    % Extract components
    G = Zi.G; A = Zi.A; b = Zi.b; c = Zi.c;

    % Tolerances (match script version)
    tol = 1e-15;
    tolNonZero = 1e-6;
    tolR = 1e-8;
    tolFeas = 1e-4; %This tolerance is checked to be absolutely needed at least for one NN run. 
    tolE = 1e-8;
    redundant_found=false;

    % Initial cleanup
    zeroCols = find(all(abs([G; A]) < tol, 1));
    G(:, zeroCols) = [];
    A(:, zeroCols) = [];

    zeroRows = find(all(abs([A, b]) < tol, 2));
    A(zeroRows, :) = [];
    b(zeroRows) = [];
    
    iterMax = 200;  % Match script version
    for iter = 1:iterMax
        [nc, ng] = size(A);
        generators_in=ng;
        constraints_in=nc;

        if nc == 0, break; end

        E = repmat([-1, 1], ng, 1);
        R = repmat([-inf, inf], ng, 1);

        for j = 1:ng

            idx = [1:j-1, j+1:ng];
            E_m = 0.5 * (E(idx,1) + E(idx,2));
            E_r = 0.5 * (E(idx,2) - E(idx,1));
            nonZeroRows = find(abs(A(:,j)) >= tolNonZero);

            for i = nonZeroRows'
                aij = A(i,j);
                middle = (1/aij) * A(i,idx) * E_m;
                range = sum(abs((1/aij) * A(i,idx) .* E_r'), 2);
                Rj_lower = b(i)/aij - (middle + range);
                Rj_upper = b(i)/aij - (middle - range);

                R(j,1) = max(R(j,1), Rj_lower);
                R(j,2) = min(R(j,2), Rj_upper);

                redundant_found=R(j,1) >= -1 + tolR && R(j,2) <= 1 - tolR;
                if redundant_found

                    %part 1 removal done here

                    [G, A, b, c] = removeGenerator(G, A, b, c, j, i);
                    break;
                else
                    E(j,1) = max(E(j,1), R(j,1));
                    E(j,2) = min(E(j,2), R(j,2));
                    if E(j,1) > E(j,2)
                        E(j,:) = E(j,[2,1]);
                    end
                end
            end
            if redundant_found, break; end
        end

        if redundant_found, continue; end

        % Rescaling
        E_m = 0.5 * (E(:,1) + E(:,2));
        E_r = 0.5 * (E(:,2) - E(:,1));
        E_r(E_r <= tolE) = 0;

        c = c + G * E_m;
        G = G * diag(E_r);
        b = b - A * E_m;
        A = A * diag(E_r);

        zeroCols = find(all(abs([G; A]) < tol, 1));
        G(:, zeroCols) = [];
        A(:, zeroCols) = [];

        zeroRows = find(all(abs([A, b]) < tol, 2));
        A(zeroRows, :) = [];
        b(zeroRows) = [];

        %Reduction only if at least zero rows or columns are found

        if isempty(zeroCols) && isempty(zeroRows)
            R(:,1) = (R(:,1) - E_m) ./ E_r;
            R(:,2) = (R(:,2) - E_m) ./ E_r;

            totalExceed = max(0, -1 - R(:,1)) + max(0, -1 + R(:,2));
            [~, indxSorted] = sort(totalExceed);

            for m = 1:ng
                j = indxSorted(m);
                redundant_found=feasibilityCheckGurobi(A, b, j, tolFeas);
                if redundant_found
                    %third part reduces if redundancy found
                    i = find(abs(A(:,j)) >= tolNonZero, 1);
                    if ~isempty(i)
                        [G, A, b, c] = removeGenerator(G, A, b, c, j, i);
                        break;
                    end
                end
            end
            if redundant_found
                continue
            end
        end
        % Additional termination condition from script
        generators_out = size(G,2);
        constraints_out=size(A,1);
        if  generators_in == generators_out && constraints_in==constraints_out
            break;
        end
    end
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

function redundant = feasibilityCheckGurobi(A, b, j, tolFeas)
    ng=size(A,2);
    s = zeros(ng,1);
    lbs = -ones(ng,1);
    ubs = ones(ng,1);
    model.sense = '=';
    model.A = sparse(A);
    model.rhs = b;
    params.Threads = 1;
    params.outputflag = 0;

    % Feasibility checks
    lbs(j) = 1+tolFeas;
    ubs(j) = inf;
    model.lb = lbs;
    model.ub = ubs;
    model.obj = s;
    result = gurobi(model,params);
    RjLPub = result.status;

    if strcmp(RjLPub,'OPTIMAL') % Would switching the order of ub vs lb reduce number of gurobi calls? 
        %yes it does matter feasibility checks for lower bounds
        %before upper bounds meant three more gurobi calls so
        % the current sequence is a more efficient one.
        RjLPlb = RjLPub;
    else
        lbs(j) = -inf;
        ubs(j) = -1-tolFeas;
        model.lb = lbs;
        model.ub = ubs;
        model.obj = -s;
        result = gurobi(model,params);
        RjLPlb = result.status;
    end
    redundant = (~strcmp(RjLPub,'OPTIMAL')) && (~strcmp(RjLPlb,'OPTIMAL'));
end