% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Checks if a given point is within a zonotope.
%       Based on Proposition 2 of Scott et al. Constrained Zonotopes: A new tool for 
%       set-based estimation and fault detection 
%   Syntax:
%       result = checkPointContain(Z, point)
%   Inputs:
%       Z - constrained zonotope
%       point - n x 1 vector defining the point to be analyzed
%       verbose (optional) - if true prints a statement about the result; default is false
%   Outputs:
%       result - returns true if point is contained within the constrained zonotope
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = checkPointContain(obj,point,varargin)

verbose = false;
if length(varargin) > 0
    verbose = true;
end

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
prob = optimproblem;
xc = optimvar('xc', obj.nG, 'LowerBound', -1, 'UpperBound', 1);

%Define objective and constraints
prob.Objective = 0; % Feasibility
prob.Constraints.Gc = obj.G*xc == point - obj.c;
prob.Constraints.Ab = obj.A*xc == obj.b;

% sol - n x 1 vector maximizing objective function subject to constraints
% fval - 1 x 1 scalar maximum value of objective function
% exitFlag - 1 x 1 scalar exit condition (see solver-specific exit condition codes)
options = optimoptions('linprog','Display','off');
[sol, fval, exitflag] = solve(prob,'Options',options);

if exitflag == 1 % solution is found so the point is contained
    out = true;
    if verbose
        disp('    RESULT: Point is contained.')
    end
else % optimization is infeasible, so the point is not contained
    out = false;
    if verbose
        disp('    RESULT: Point is not contained.')
    end
end

end