% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the value of the support function for a given set and
%       direction
%   Syntax:
%       [s,x] = supportFunc(X,dims,d)
%   Inputs:
%       X - zonotopic set in R^n (memZono object)
%       dims - dimLabels for applying support vector to
%       d 
%   Outputs:
%       s - scalar such that s = max(d'*x), where x \in X
%       x - n x 1 vector such that x = argmax(d'*x), where x \in X
%   Notes:
%       Defined to find maximum in desired dimension
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [s,x] = supportFunc(obj,dims,d)
    arguments
        obj memZono
        dims
        d = [];
    end

    if isempty(d)
        d = ones(length(dims),1);
    end

    if isnumeric(obj)
        [s,x] = obj(dims).Z.supportFunc(d);
    else
        [s,x] = supportFunSpecial(obj(dims),d);
    end

end


function [s,x] = supportFunSpecial(obj,d)
    

end