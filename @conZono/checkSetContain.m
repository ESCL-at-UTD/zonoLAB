
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Checks if the constrained zonotope X is contained within the constrained zonotope Y.
%   Syntax:
%       [sol,fval,exitFlag] = checkSetContain(X, Y, verbose)
%   Inputs:
%       X - constrained zonotope
%       Y - constrained zonotope  
%       verbose (optional) - if true prints a statement about the result; default is false
%   Outputs:
%       result - returns true if X is contained within Y, false otherwise.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [out] = checkSetContain(X,Y,varargin)

verbose = false;
if length(varargin) > 0
    verbose = true;
end

%Convert ConZonos to AH Polytopes
[X_AH.G, X_AH.c, X_AH.H, X_AH.k] = conZono2AHPoly(X);
[Y_AH.G, Y_AH.c, Y_AH.H, Y_AH.k] = conZono2AHPoly(Y);

n_x = size(X_AH.G, 2);
n_y = size(Y_AH.G, 2);
n_hx = size(X_AH.H, 1);
n_hy = size(Y_AH.H, 1);

% Create optimization variables
prob = optimproblem;
gamma = optimvar('gamma', n_y, n_x);
beta = optimvar('beta', n_y);
lambda = optimvar('lambda', n_hy, n_hx);
l = optimvar('l', n_hy, n_hx);

prob.Objective = 0; %Feasibility
prob.Constraints.gamma = X_AH.G == Y_AH.G*gamma;
prob.Constraints.beta = Y_AH.c - X_AH.c == Y_AH.G*beta;
prob.Constraints.lambda_abs_1 = l >= lambda;
prob.Constraints.lambda_abs_2 = l >= -lambda;
prob.Constraints.lambda1 = l*X_AH.H == Y_AH.H*gamma;
prob.Constraints.lambda2 = l*X_AH.k <= Y_AH.k + Y_AH.H*beta;

options = optimoptions('linprog','Display','off');
[sol, fval, exitflag] = solve(prob, 'Options',options);

if exitflag == 1 % solution is found, so the set is contained
    out = true;
    if verbose
        disp('    RESULT: Set is contained.')
    end
else % optimization is infeasible, so the set is not contained
    out = false;
    if verbose
        disp('    RESULT: Set is not contained.')
    end
end

end