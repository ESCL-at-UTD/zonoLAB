
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

%create problem to pass to enforceSetContain
prob = optimproblem;
prob.Objective = 0; %Feasibility
prob = enforceSetContain(X,Y,prob);

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