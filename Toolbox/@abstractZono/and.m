% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the generalized intersection of two sets, 
%       Z = X \cap_R Y = {x \in X | Rx \in Y}
%   Syntax:
%       Z = and(X,Y,R)        for generalized intersection
%       Z = and(X,Y) = X & Y  for regular intersection (R = I)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       Y - zonotopic set in R^m (hybZono, conZono, or zono object)
%       R - m x n real matrix
%   Outputs:
%       Z - zonotopic set in R^n (hybZono or conZono object)
%   Notes:
%       Overloaded '&' operator for regular intersection (M = I)
%       If X or Y is hybZono, Z is hybZono
%       If X and Y zono, Z is conZono
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = and(obj1,obj2,R)

if nargin == 2
    R = eye(obj1.n);
end

if (obj1.n ~= size(R,2)) || (obj2.n ~= size(R,1))
    error('Inconsistent dimensions.')
end

% Ensure object classes match
[obj1, obj2] = matchSetType(obj1,obj2);

c = obj1.c;
switch class(obj1)
    case 'zono'
        G = [obj1.G zeros(obj1.n,obj2.nG)];
        A = [R*obj1.G -obj2.G];
        b = [obj2.c-R*obj1.c];
        out = conZono(G,c,A,b);
    case 'conZono'
        G = [obj1.G zeros(obj1.n,obj2.nG)];
        A = blkdiag(obj1.A,obj2.A);
        A = [A; R*obj1.G -obj2.G];
        b = [obj1.b; obj2.b; obj2.c-R*obj1.c];
        out = conZono(G,c,A,b);
    case 'hybZono'
        Gc = [obj1.Gc zeros(obj1.n,obj2.nGc)];
        Gb = [obj1.Gb zeros(obj1.n,obj2.nGb)];
        Ac = blkdiag(obj1.Ac,obj2.Ac);
        Ac = [Ac; R*obj1.Gc -obj2.Gc];
        Ab = blkdiag(obj1.Ab,obj2.Ab);
        Ab = [Ab; R*obj1.Gb -obj2.Gb];
        b = [obj1.b; obj2.b; obj2.c-R*obj1.c];
        out = hybZono(Gc,Gb,c,Ac,Ab,b);
end

end