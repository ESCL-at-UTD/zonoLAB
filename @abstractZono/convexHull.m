% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the convex hull of two zonotopic sets, Z = CH(X \cup Y), or
%           the convex hull of a single zonotopic set. 
%   Syntax:
%       Z = convexHull(X,Y)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       Y - zonotopic set in R^n (hybZono, conZono, or zono object)
%   Outputs:
%       Z - zonotopic set in R^n (conZono object)
%   Notes:
%       Current method for convex hull of hybZono introduces exponentially
%       many new constraints and continuous variables; see the function
%       sharpHybZono for more detail. 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = convexHull(obj1,varargin)

switch nargin
    case 2
        obj2 = varargin{1};
    case 1
        if isa(obj1, 'zono')
            out = conZono(obj1);
        elseif isa(obj1, 'conZono')
            out = obj1;
        elseif isa(obj1, 'hybZono')
            out = convexRelaxation(sharpHybZono(obj1));
        end
        return
    otherwise
        error('The method convexHull accepts only one or two inputs.')
end

        if isa(obj1,'zono')
            obj1 = conZono(obj1);
        end
        if isa(obj2,'zono')
            obj2 = conZono(obj2);
        end
        if isa(obj1, 'hybZono')
            obj1 = convexRelaxation(sharpHybZono(obj1));
        end
        if isa(obj2, 'hybZono')
            obj2 = convexRelaxation(sharpHybZono(obj2));
        end

        c = (obj1.c + obj2.c)/2;
        G = [obj1.G obj2.G (obj1.c-obj2.c)/2 zeros(obj1.n,2*(obj1.nG+obj2.nG))];

        H = [eye(obj1.nG) zeros(obj1.nG,obj2.nG)  -0.5*ones(obj1.nG,1);...
            -eye(obj1.nG) zeros(obj1.nG,obj2.nG)  -0.5*ones(obj1.nG,1);...
             zeros(obj2.nG,obj1.nG)  eye(obj2.nG)  0.5*ones(obj2.nG,1);...
             zeros(obj2.nG,obj1.nG) -eye(obj2.nG)  0.5*ones(obj2.nG,1)];
        I = eye(2*(obj1.nG + obj2.nG));
        f = -0.5*ones(2*(obj1.nG + obj2.nG),1);

        A = [obj1.A zeros(obj1.nC,obj2.nG) -obj1.b/2 zeros(obj1.nC,size(G,2)-obj1.nG-obj2.nG-1);...
            zeros(obj2.nC,obj1.nG) obj2.A   obj2.b/2 zeros(obj2.nC,size(G,2)-obj2.nG-obj1.nG-1);...
            H I];
        b = [obj1.b/2;obj2.b/2;f];

        out = conZono(G,c,A,b);
end