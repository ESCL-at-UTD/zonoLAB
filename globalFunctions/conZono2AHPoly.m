% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Transforms a constrained zonotope to an AH-polytope (affine
%       H-polytope) 
%   Syntax:
%       result = conZono2AHPoly(Z)
%   Inputs:
%       Z - constrained zonotope
%   Outputs:
%       [G, c, H, k] - returns the AH-polytope representation matrices of Z.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [G,c,H,k] = conZono2AHPoly(Z)

if Z.nC == 0
    G = Z.G;
    c = Z.c;    
    H = [-eye(Z.nG);eye(Z.nG)];
    k = ones(2*Z.nG,1);
else
    T = null(Z.A,'r');  % Computes the right null space of A.
    s = pinv(Z.A)*Z.b;
    P = Polyhedron('H',[T ones(size(s,1),1)-s; -T ones(size(s,1),1)+s]); % Inequality constraints on \xi variables.
    P = minHRep(P); % Removes redundant halfspaces
    G = Z.G*T;
    c = Z.c+Z.G*s;    
    H = P.H(:,1:end-1);
    k = P.H(:,end);
end

end