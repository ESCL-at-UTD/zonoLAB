% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the axis-aligned bounding-box of a set
%   Syntax:
%       [B] = boundingBox(X)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%   Outputs:
%       B - zonotope in R^n with diagonal generator matrix
%   Notes:
%        X \subseteq B
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [out] = boundingBox(obj)

lb = zeros(obj.n,1);
ub = zeros(obj.n,1);
basisVectors = eye(obj.n);
for i = 1:obj.n
    [s,~] = supportFunc(obj,-basisVectors(:,i));
    lb(i) = -s;
    [s,~] = supportFunc(obj,basisVectors(:,i));
    ub(i) = s;
end

out = zono(diag(ub-lb)/2,(ub+lb)/2);

end