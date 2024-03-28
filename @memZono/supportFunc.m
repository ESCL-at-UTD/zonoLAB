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
        [s,x] = projection(obj,dims).Z.supportFunc(d);
    else
        if issym(obj)
            
        else %<=== assume optimvar
            [s,x] = supportFunOptimvar(projection(obj,dims),d);
        end
    end

end


function [s,x] = supportFunOptimvar(obj,d)
    switch obj.baseClass
        case 'zono'
            xi = fcn2optimexpr(@(G) sign(d'*G)', obj.G);
            x = obj.c + obj.G*xi;
            s = d'*x;

        case 'conZono'
            ub = ones(obj.nG,1); lb = -ub;
            xi = fcn2optimexpr(@(G,A,b) ...
                linprog(-(d'*G)',[],[],A,b,lb,ub),...
                obj.G,obj.A,obj.b);
            x = obj.c + obj.G*xi;
            s = d'*x;

        case 'hybZono'
            error('supportFun not currently implimented for optimvar')
        
    end
end