function [sol, fval, exitflag] = checkPointContain(point,zono,verbose)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Global function:
%       Check if a point is within a zonotope.
%       Based on Proposition 2 of Scott et al. Constrained Zonotopes: A new tool for 
%       set-based estimation and fault detection 
%   Syntax:
%       [sol,fval,exitFlag] = pointContain(point, zono, verbose)
%   Inputs:
%       point - n x 1 vector defining point to be analyzed
%       zono - zonotope, constrained zonotope, or hybrid zonotope  
%       verbose - true to send message to console if point is contained or not,
%                 false for no message confirmation  
%   Outputs:
%       sol - n x 1 vector maximizing objective function subject to constraints
%       fval - 1 x 1 scalar maximum value of objective function
%       exitFlag - 1 x 1 scalar exit condition (see solver-specific exit condition codes)
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %


%Create optimization variables
pointContain = optimproblem;
xc = optimvar('xc', zono.nG, 'LowerBound', -1, 'UpperBound', 1);

obj = 0; %Feasibility

%Check for constraint matrices
%Have to create an instance of A and b in the function if it does not exist
%If input is a regular zonotope, set constraints to zero
A = zono.A;
b = zono.b;
if isempty(zono.A)
    zono.A = 0;
end
if isempty(zono.b)
    zono.b = 0;
end

%Define objective and constraints
pointContain.Objective = obj;
pointContain.Constraints.Gc = zono.G*xc == point - zono.c;
pointContain.Constraints.Ab = zono.A*xc == zono.b;

[sol, fval, exitflag] = solve(pointContain);

if verbose == true
    if exitflag == 1
        fprintf("The point is contained.\n");
    else
        fprintf("The point is not contained.\n");
    end
end
end