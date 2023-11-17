% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the simplest zonotopic set representation
%   Syntax:
%       [Z] = simplifySetType(X)
%   Inputs:
%       X - zonotopic set (hybZono, conZono, or zono object)
%   Outputs:
%       Z - zonotopic set (hybZono, conZono, or zono object)
%   Notes:
%       A conZono without constraints is a zono
%       A hybZono without binary factors and without constraints is a zono
%       A hybZono without binary factors but with constraints is a conZono
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [out] = simplifySetType(obj)

% Keep the original object unless:
out = obj;

if isa(obj,'conZono') && obj.nC == 0 % a conZono without constraints is a zono
    out = zono(obj);
elseif isa(obj,'hybZono') && obj.nGb == 0 % a hybZono without binary factors
    if obj.nC == 0 % and without constraints is a zono
        out = zono(obj);
    else % otherwise a conZono
        out = conZono(obj);
    end
end

end