% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns two zonotopic sets as the same type (hybZono, conZono, or zono object)
%   Syntax:
%       [Z1, Z2] = matchSetType(X1,X2)
%   Inputs:
%       X1 - zonotopic set (hybZono, conZono, or zono object or vector)
%       X2 - zonotopic set (hybZono, conZono, or zono object or vector)
%   Outputs:
%       Z1 - zonotopic set (hybZono, conZono, or zono object)
%       Z2 - zonotopic set (hybZono, conZono, or zono object)
%   Notes:
%       Z1 and Z2 will always be the same class
%       If X1 or X2 is hybZono, then Z1 and Z2 are hybZono
%       If X1 or X2 is conZono, then Z1 and Z2 are conZono
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [out1, out2] = matchSetType(obj1,obj2)

if strcmp(class(obj1),class(obj2)) % If classes match
    out1 = obj1;
    out2 = obj2;
    return
end

% Continue if classes do not match

types = [isa(obj1,'hybZono') isa(obj2,'hybZono');...
         isa(obj1,'conZono') isa(obj2,'conZono');...
         isa(obj1,'zono')    isa(obj2,'zono')];
if sum(types(1,:)) ~=0 % One of the sets is a hybrid zonotope
    if types(1,1) == 1
        out1 = obj1;
        out2 = hybZono(obj2);
    else
        out1 = hybZono(obj1);
        out2 = obj2;
    end
elseif sum(types(2,:)) ~=0 % One of the sets is a constrained zonotope
    if types(2,1) == 1
        out1 = obj1;
        out2 = conZono(obj2);
    else
        out1 = conZono(obj1);
        out2 = obj2;
    end
elseif sum(types(3,:)) ~=0 % One of the sets is a zonotope
    if types(3,1) == 1
        out1 = obj1;
        out2 = zono(obj2);
    else
        out1 = zono(obj1);
        out2 = obj2;
    end
else
    warning(['At least on set must be a hybrid zonotopes,';...
        'constrained zonotopes, or zonotopes.'])
end

end