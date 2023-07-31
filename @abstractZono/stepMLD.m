% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the 1-step reachable set for an MLD system.
%   Syntax:
%       Z = stepMLD(X)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       
%   Outputs:
%       Z - hybZono in R^m
%   Notes:
%       MLD of the form:
%           x^+ =   A x + B_u u + B_w w +  B_aff
%           s.t.  E_x x + E_u u + E_w w <= E_aff
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = stepMLD(X0,U,W,A,B_u,B_w,B_aff,E_x,E_u,E_w,E_aff)

nE = size(E_aff,1);     % Number of inequality constraints for MLD

if isempty(U)
    U = zono(1);
    B_u = zeros(X0.n,1);
    E_u = zeros(nE,1);
end

V = [B_u; E_u]*U + [B_w; E_w]*W + [B_aff; zeros(nE,1)];

Y = [A; E_x]*X0 + V;

H = eye(nE);
f = E_aff;
R = [zeros(nE,X0.n) eye(nE)];

out = [eye(X0.n) zeros(X0.n,nE)]*halfspaceIntersection(Y,H,f,R);

end