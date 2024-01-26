% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Checks if the zonotope X is contained within a zonotope Y.
%   Syntax:
%       [sol,fval,exitFlag] = setContain(X, Y, verbose)
%   Inputs:
%       X - constrained zonotope
%       Y - constrained zonotope  
%       verbose - true to send message to console if point is contained or not,
%                 false for no message confirmation  
%   Outputs:
%       sol - n x 1 vector maximizing objective function subject to constraints
%       fval - 1 x 1 scalar maximum value of objective function
%       exitFlag - 1 x 1 scalar exit condition (see solver-specific exit condition codes)
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [sol, fval, exitflag] = setContain(X, Y, verbose)

setContain = optimproblem;

obj = 0; %Feasibility

n_x = size(X.G, 2);
n_y = size(Y.G, 2);

%Define Optimization variables
gamma = optimvar('gamma', n_y, n_x);
beta = optimvar('beta', n_y);

%Initialize variables for nonlinear function calculation
initial.gamma = zeros(n_y, n_x);
initial.beta = zeros(n_y,1);

setContain.Objective = obj;
setContain.Constraints.gamma = X.G == Y.G*gamma;
setContain.Constraints.beta = Y.c - X.c == Y.G*beta;
setContain.Constraints.gamma_beta = fcn2optimexpr(@abs, gamma)*ones(n_x,1) + fcn2optimexpr(@abs, beta) <= ones(n_y,1); 

% sol - n x 1 vector maximizing objective function subject to constraints
% fval - 1 x 1 scalar maximum value of objective function
% exitFlag - 1 x 1 scalar exit condition (see solver-specific exit condition codes)
[sol, fval, exitflag] = solve(setContain,initial);
if verbose == true
    if exitflag == 1
        fprintf("The set is contained.\n");
    else
        fprintf("The set is not contained.\n");
    end
end

end