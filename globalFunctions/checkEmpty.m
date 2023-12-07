function [sol, fval, exitflag] = checkEmpty(zono,verbose)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Global function:
%       Check if a zonotope is an empty set.
%       Based on Proposition 2 of Scott et al. Constrained Zonotopes: A new tool for 
%       set-based estimation and fault detection 
%   Syntax:
%       [sol,fval,exitFlag] = zonoEmpty(zono, verbose)
%   Inputs:
%       zono - constrained zonotope or hybrid zonotope  
%       verbose - true to send message to console if zonotope is empty or not,
%                 false for no message confirmation  
%   Outputs:
%       sol - n x 1 vector maximizing objective function subject to constraints
%       fval - 1 x 1 scalar maximum value of objective function
%       exitFlag - 1 x 1 scalar exit condition (see solver-specific exit condition codes)
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%Create optimization variables
isEmpty = optimproblem;
xc = optimvar('xc', zono.nG, 'LowerBound', -1, 'UpperBound', 1);

obj = 0; %Feasibility

%Check for constraint matrices
%Have to create an instance of A and b in the function if it does not exist
A = zono.A;
b = zono.b;
%If input is a regular zonotope, set constraints to zero
if isempty(zono.A)
    zono.A = 0;
end
if isempty(zono.b)
    zono.b = 0;
end

%Define objective and contstraints
isEmpty.Objective = obj;
isEmpty.Constraints.Ab = zono.A*xc == zono.b;

[sol, fval, exitflag] = solve(isEmpty);

if verbose == true
    if exitflag == 1
        fprintf("The set is not empty.\n");
    else
        fprintf("The set is empty.\n");
    end
end

end