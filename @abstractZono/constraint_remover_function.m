function [G, c, A, b] = constraint_remover_function(G, c, A, b, remove_position,removable_constraint)
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
    E_lam = zeros(n_g, n_c);
    E_lam(remove_position, removable_constraint) = 1; % Place a 1 at the (remove_position, removable_constraint) position

    % Step 2: Compute \Lambda_G and \Lambda_A
    a_1j = A(removable_constraint, remove_position); % The removable constraint row
    Lambda_G = G * E_lam / a_1j; % \Lambda_G = G * E * a_1j^(-1)
    Lambda_A = A * E_lam / a_1j; % \Lambda_A = A * E * a_1j^(-1)

    % Step 3: Update the generator matrix G
    G = G - Lambda_G * A;

    % Step 4: Update the center vector c
    c = c + Lambda_G * b;

    % Step 5: Update the constraint matrix A
    A = A - Lambda_A * A;

    % Step 6: Update the constraint vector b
    b = b - Lambda_A * b;

    % Step 7: Remove the selected constraint from all matrices
    G(:, remove_position) = []; % Remove the column in G corresponding to remove_position
    A(removable_constraint, :) = []; % Remove the row in A corresponding to removable_constraint
    A(:, remove_position) = []; % Remove the column in A corresponding to remove_position
    b(removable_constraint,:) = [];
    end

