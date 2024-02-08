% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the Pontryagin difference of two zonotopic sets, Z = X \ominus Y
%   Syntax:
%       Z = pontryDiff(X,Y)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       Y - zonotope in R^n (zono object)
%   Outputs:
%       Z - zonotopic set in R^n (hybZono or conZono object)
%   Notes:
%       For exact computation, Y must be a zono object.
%       Approximate methods are not implemented yet.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = pontryDiff(obj1,obj2)

out = obj1;
out.c = obj1.c - obj2.c;

for i = 1:obj2.nG
    out_plus = out;
    out_plus.c = out.c + obj2.G(:,i);
    out_minus = out;
    out_minus.c = out.c - obj2.G(:,i);
    out = out_plus & out_minus;
end