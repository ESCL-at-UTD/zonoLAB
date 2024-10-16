% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the n-1 dimensional constrained zonotope representation of 
%       an axis-aligned cross-section of an n-dimensional zono or conZono object.
%   Syntax:
%       Z = slice(obj,R,v)
%   Inputs:
%       obj - zonotopic set in R^n (conZono or zono object; hybZono not currently supported)
%       R - Row vector (1 x R^n) representing the normal vector of the
%       slicing hyperplane.
%       v - Scalar offset of the slicing hyperplane
%   Outputs:
%       Z - conZono object in R^n 
%       (Slice operation adds one continuous constraint, so the slice of a 
%       zono is a conZono, and the slice of a conZono is a conZono 
%       with one more constraint Rx = v.)
%   Notes:
%       
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = zonoSlice(obj,R,v)

G = obj.G; c = obj.c;

switch class(obj)
    case 'zono'
        A = R*obj.G;
        b = v-(R*obj.c);
        out = conZono(G,c,A,b);
    case 'conZono'
        A = [obj.A; R*obj.G];
        b = [obj.b; v-(R*obj.c)];
        out = conZono(G,c,A,b);
    case 'hybZono'
        Gb = obj.Gb;
        Ac = [obj.Ac; R*obj.Gc];
        Ab = [obj.Ab; R*obj.Gb];
        b  = [obj.b; v-(R*obj.c)];
        out = hybZono(G,Gb,c,Ac,Ab,b);
end



% R = zeros(1, obj.n); R(dim) = 1;
% switch class(obj)
%     case 'zono'
%         A = R*G; b = v-(R*c);
%     case 'conZono'
%         A = [obj.A; R*G]; b = [obj.b; v-(R*c)];
%     case 'hybZono'
%         error('Hybrid Zonotope input objects are not currently supported with this function.')
% end
% 
% [~,~,exit_flag] = linprog(zeros(1,obj.nG),[],[],A,b, -1*ones(obj.nG,1), ones(obj.nG,1));
% if exit_flag == 1
%     Y = conZono(G,c,A,b);
%     projection_dims = 1:obj.n; 
%     projection_dims(dim) = [];
% else
%     error('The output conZono object will be empty, because the cross-sectional slice described by inputs dim and v does not intersect the input zono/conZono obj.');
% end
% out = projection(Y,projection_dims);

end