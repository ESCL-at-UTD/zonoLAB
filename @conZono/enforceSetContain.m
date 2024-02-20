% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Enforces constraints to check if the constrained zonotope X is contained within the constrained zonotope Y.
%   Syntax:
%       [prob] = enforceSetContain(X, Y, prob, varargin)
%   Inputs:
%       X - constrained zonotope
%       Y - constrained zonotope
%       prob - MATLAB optimproblem
%       varargin - string or character input to add as a suffix for variables if doing nested optimization to
%       enforce uniqueness (default is '1'.)
%
%       
%   Outputs:
%       prob - returns the optimproblem containing the constraints necessary to check for constrained zonotope set containment.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [prob] = enforceSetContain(X,Y,prob,varargin)

%check for user input of suffix variable, otherwise 1
if length(varargin) > 0
    suffix = string(varargin);
else
    suffix = '1';
end

%Convert ConZonos to AH Polytopes
[X_AH.G, X_AH.c, X_AH.H, X_AH.k] = conZono2AHPoly(X);
[Y_AH.G, Y_AH.c, Y_AH.H, Y_AH.k] = conZono2AHPoly(Y);

n_x = size(X_AH.G, 2);
n_y = size(Y_AH.G, 2);
n_hx = size(X_AH.H, 1);
n_hy = size(Y_AH.H, 1);

% Create optimization variables
gamma = optimvar(sprintf('gamma_%s',suffix), n_y, n_x);
beta = optimvar(sprintf('beta_%s',suffix), n_y);
lambda = optimvar(sprintf('lambda_%s',suffix), n_hy, n_hx);

%if constraints structure is empty, create structure with first constraint
if isempty(prob.Constraints)
    prob.Constraints = struct(sprintf('gamma_%s',suffix),X_AH.G == Y_AH.G*gamma);
else
    prob.Constraints = setfield(prob.Constraints, sprintf('gamma_%s',suffix), X_AH.G == Y_AH.G*gamma) ;
end

prob.Constraints = setfield(prob.Constraints, sprintf('beta_%s',suffix), Y_AH.c - X_AH.c == Y_AH.G*beta) ;
prob.Constraints = setfield(prob.Constraints, sprintf('lambda_nonneg_%s',suffix), lambda >= 0) ;
prob.Constraints = setfield(prob.Constraints, sprintf('lambda_gamma_%s',suffix), lambda*X_AH.H == Y_AH.H*gamma) ;
prob.Constraints = setfield(prob.Constraints, sprintf('lambda_beta_%s',suffix), lambda*X_AH.k <= Y_AH.k + Y_AH.H*beta);

end

