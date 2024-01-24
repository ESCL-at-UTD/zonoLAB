% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Checks if a zonotope is an empty set.
%       Based on Proposition 2 of Scott et al. Constrained Zonotopes: A new tool for 
%       set-based estimation and fault detection 
%   Syntax:
%       result = checkEmpty(Z)
%   Inputs:
%       Z - constrained zonotope
%   Outputs:
%       result - returns true if the set is empty, otherwise false
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = checkEmpty(obj)

% Check for constraint matrices
% Have to create an instance of A and b in the function if it does not exist
% If input is a regular zonotope, set constraints to zero
if isempty(obj.A)
    obj.A = 0;
end
if isempty(obj.b)
    obj.b = 0;
end

% Create optimization variables
model = optimproblem;
xc = optimvar('xc', obj.nG, 'LowerBound', -1, 'UpperBound', 1);

% Define objective and contstraints
model.Objective = 0; % Feasibility
model.Constraints.Ab = obj.A*xc == obj.b;

% sol - n x 1 vector maximizing objective function subject to constraints
% fval - 1 x 1 scalar maximum value of objective function
% exitFlag - 1 x 1 scalar exit condition (see solver-specific exit condition codes)
[sol, fval, exitflag] = solve(model);

if exitflag == 1 % optimization is feasible and so the set is not empty
    out = false;
else % optimization is infeasible and so the set is empty
    out = true;
end

end