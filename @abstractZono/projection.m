% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the projection of a set onto the specified dimensions.
%   Syntax:
%       Z = projection(X,dims)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       dims - 1 x m vector (m <= n) with the desired dimensions
%   Outputs:
%       Z - zonotopic set in R^m (hybZono, conZono, or zono object)
%   Notes:
%       If X \in R^3, dims = [1 3] will project X onto the x and z axes.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = projection(obj,dims)

out = obj;

out.c = obj.c(dims,:);
switch class(obj)
    case 'zono'
        out.G = obj.G(dims,:);
    case 'conZono'
        out.G = obj.G(dims,:);
    case 'hybZono'
        out.Gc = obj.Gc(dims,:);
        out.Gb = obj.Gb(dims,:);
end

end