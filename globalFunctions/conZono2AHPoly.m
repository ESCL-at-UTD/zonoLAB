% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Transforms a constrained zonotope to an AH-polytope (affine
%       H-polytope) 
%   Syntax:
%       result = conZono2AHPoly(Z)
%   Inputs:
%       Z - constrained zonotope
%   Outputs:
%       result - returns the AH-polytope representation of Z.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [result] = conZono2AHPoly(Z)

if Z.nC == 0
    result.G = Z.G;
    result.c = Z.c;    
    result.A = [-eye(Z.nG);eye(Z.nG)];
    result.b = ones(2*Z.nG,1);
else
    T = null(Z.A,'r');  % Computes the right null space of A.
    s = pinv(Z.A)*obj.b;
    P = Polyhedron('H',[T ones(size(s,1),1)-s; -T ones(size(s,1),1)+s]); % Inequality constraints on \xi variables.
    P = minHRep(P); % Removes redundant halfspaces
    result.G = Z.G*T;
    result.c = Z.c+Z.G*s;    
    result.A = P.H(:,1:end-1);
    result.b = P.H(:,end);
end

end