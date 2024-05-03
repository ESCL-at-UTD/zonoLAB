% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the intersection of zonotopic set X with halfspace(s) H Rx <= f
%   Syntax:
%       Z = halfspaceIntersection(X,H,f,R) for generalized intersection
%       Z = halfspaceIntersection(X,H,f)   for regular intersection (R = I)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       H - nH x m real matrix defining the nH halfspace normal vectors
%       f - nH x 1 real vector defining the nH offsets
%       R - m x n real matrix
%   Outputs:
%       Z - zonotopic set in R^n (hybZono, conZono, or zono object)
%   Notes:
%       If Z is zono, then X is conZono (otherwise Z and X are the same
%       type).
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = halfspaceIntersection(obj,H,f,R)

if nargin <= 3
    R = eye(obj.n);
end

% Check for compatable dimensions
if size(H,2) ~= size(R,1)
    error(['Second argument must be a vector or matrix with ',num2str(size(R,1)),' columns.'])
elseif size(R,2) ~= obj.n
    error(['Fourth argument must be a vector or matrix with ',num2str(obj.n),' columns.'])
elseif size(H,1) ~= size(f,1)
    error('Second and third arguments must have the same number of rows.')
elseif size(f,2) ~= 1
    error('Third argument must be a column vector.')
end

nH = size(H,1);
switch class(obj)
    case 'zono'
        obj = conZono(obj);
        out = halfspaceIntersection(obj,H,f,R);
    case 'conZono'
        G = obj.G;
        c = obj.c;
        A = obj.A;
        b = obj.b;
        d_max = f-H*R*c + sum(abs(H*R*G),2); % Computes the distance of the hyperplane from the farthest vertex of the zonotope
        A = [A zeros(size(A,1),nH); H*R*G diag(d_max)/2];
        b = [b; f-H*R*c - d_max/2];
        G = [G zeros(obj.n,nH)];
        out = conZono(G,c,A,b);
    case 'hybZono'
        Gc = obj.Gc;
        Gb = obj.Gb;
        c = obj.c;
        Ac = obj.Ac;
        Ab = obj.Ab;
        b = obj.b;
        d_max = f-H*R*c + sum(abs(H*R*[Gc Gb]),2); % Computes the distance of the hyperplane from the farthest vertex of the zonotope
        Ac = [Ac zeros(size(Ac,1),nH); H*R*Gc diag(d_max)/2];
        Ab = [Ab; H*R*Gb];
        b = [b; f-H*R*c - d_max/2];
        Gc = [Gc zeros(obj.n,nH)];
        out = hybZono(Gc,Gb,c,Ac,Ab,b);
end

end