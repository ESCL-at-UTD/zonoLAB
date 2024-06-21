% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Checks if the zonotope X is contained within the zonotope Y.
%   Syntax:
%       [sol,fval,exitFlag] = checkSetContain(X, Y, verbose)
%   Inputs:
%       X - zonotope
%       Y - zonotope  
%       verbose (optional) - if true prints a statement about the result; default is false
%   Outputs:
%       result - returns true if X is contained within Y, false otherwise.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [out] = checkSetContain(obj, Y, varargin)
X = obj;

verbose = false;
if length(varargin) > 0
    verbose = true;
end

n_x = size(X.G, 2);
n_y = size(Y.G, 2);

% Create optimization variables
prob = optimproblem;
gamma = optimvar('gamma', n_y, n_x);
beta = optimvar('beta', n_y);
g = optimvar('g', n_y, n_x);
b = optimvar('b', n_y);

% Initialize variables for nonlinear function calculation
%initial.gamma = ones(n_y, n_x)/(n_x+2);
%initial.beta = ones(n_y,1)/(n_x+2);

prob.Objective = 0; % Feasibility
prob.Constraints.gamma = X.G == Y.G*gamma;
prob.Constraints.beta = Y.c - X.c == Y.G*beta;
%prob.Constraints.gamma_beta = fcn2optimexpr(@abs, gamma)*ones(n_x,1) + fcn2optimexpr(@abs, beta) <= ones(n_y,1); 
%prob.Constraints.gamma_beta = sqrt(gamma.^2)*ones(n_x,1) + sqrt(beta.^2) <= ones(n_y,1); 
prob.Constraints.gamma_abs_1 = g >= gamma;
prob.Constraints.gamma_abs_2 = g >= -gamma;
prob.Constraints.beta_abs_1 = b >= beta;
prob.Constraints.beta_abs_2 = b >= -beta;
prob.Constraints.gamma_beta = g*ones(n_x,1) + b <= ones(n_y,1);

% sol - n x 1 vector maximizing objective function subject to constraints
% fval - 1 x 1 scalar maximum value of objective function
% exitFlag - 1 x 1 scalar exit condition (see solver-specific exit condition codes)
options = optimoptions('linprog','Display','off');
% [sol, fval, exitflag] = solve(prob,initial,'Options',options);
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