% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the Cartesian product of two zonotopic sets, Z = X x Y
%   Syntax:
%       Z = cartProd(X,Y)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       Y - zonotopic set in R^m (hybZono, conZono, or zono object)
%   Outputs:
%       Z - zonotopic set in R^(n+m) (hybZono, conZono, or zono object)
%   Notes:
%       If X is a hybZono and Y is a conZono, then Z is a hybZono
%       If X is a conZono and Y is a zono, then Z is a conZono
%       If X and Y are zono, then Z is a zono
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = cartProd(obj1,obj2)

% Ensure object classes match
[obj1, obj2] = matchSetType(obj1,obj2);

switch class(obj1)
    case 'zono'
        G = blkdiag(obj1.G,obj2.G);
        c = [obj1.c; obj2.c];
        out = zono(G,c);
    case 'conZono'
        G = blkdiag(obj1.G,obj2.G);
        c = [obj1.c; obj2.c];
        A = blkdiag(obj1.A,obj2.A);
        b = [obj1.b; obj2.b];
        out = conZono(G,c,A,b);
    case 'hybZono'
        Gc = blkdiag(obj1.Gc,obj2.Gc);
        Gb = blkdiag(obj1.Gb,obj2.Gb);
        c = [obj1.c; obj2.c];
        Ac = blkdiag(obj1.Ac,obj2.Ac);
        Ab = blkdiag(obj1.Ab,obj2.Ab);
        b = [obj1.b; obj2.b];
        out = hybZono(Gc,Gb,c,Ac,Ab,b);
end

end