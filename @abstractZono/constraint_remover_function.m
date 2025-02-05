function [G_new, c_new, A_new, b_new] = constraint_remover_function(G, c, A, b, remove_position,removable_constraint)
% CONSTRAINT_REDUCTION - Perform constraint reduction for zonotopes.
% 
% Inputs:
%   G - Generator matrix (n x ng)
%   c - Center vector (n x 1)
%   A - Constraint matrix (n_c x n_g)
%   b - Constraint vector (n_c x 1)
%   removable_constraint - Index of the constraint to be removed (scalar)
%
% Outputs:
%   G_new - Updated generator matrix (n x ng)
%   c_new - Updated center vector (n x 1)
%   A_new - Updated constraint matrix (n_c-1 x n_g)
%   b_new - Updated constraint vector (n_c-1 x 1)

    % Validate inputs
    [n_c, n_g] = size(A);

    % Step 1: Define the E matrix
    E = zeros(n_g, n_c);
    E(remove_position, removable_constraint) = 1; % Place a 1 at the (j, 1) position

    % Step 2: Compute \Lambda_G and \Lambda_A
    A_j = A(removable_constraint, :); % The j-th row of A
    a_1j = A_j(remove_position); % The (j, 1) element of A
    Lambda_G = G * E / a_1j; % \Lambda_G = G * E * a_1j^(-1)
    Lambda_A = A * E / a_1j; % \Lambda_A = A * E * a_1j^(-1)

    % Step 3: Update the generator matrix G
    G_new = G - Lambda_G * A;
    
    % Step 4: Update the center vector c
    c_new = c + Lambda_G * b;

    % Step 5: Update the constraint matrix A
    A_new = A - Lambda_A * A;

    % Step 6: Update the constraint vector b
    b_new = b - Lambda_A * b;

    G_new(:,remove_position)=[];
    A_new(removable_constraint, :) = [];
    A_new(:,remove_position)=[];
    b_new(removable_constraint) = [];

    % Display results (optional for debugging)
    % disp('Updated Generator Matrix G:');
    % disp(G_new);
    % disp('Updated Center Vector c:');
    % disp(c_new);
    % disp('Updated Constraint Matrix A:');
    % disp(A_new);
    % disp('Updated Constraint Vector b:');
    % disp(b_new);
    end

