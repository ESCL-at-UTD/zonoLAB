% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the complement of a set within specified bounds.
%   Syntax:
%       Z = complement(X,bounds)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       bounds - n x 1 vector bounds vector so that X \subset [-bounds,bounds]
%   Outputs:
%       Zc - zonotopic set in R^n (hybZono object)
%   Notes:
%       Z is technically the closure of the complement of X within the 
%       interval region of interest [-bounds,bounds] 
%       so that Zc = [-bounds,bounds] \ X
%       If X is a hybZono object, then it is decomposed into conZono objects and
%       the complement is returned as the intersection over the complements of
%       all nonempty conZono objects
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = complement(obj,bounds)

optSolver = solverOptions; % Using default options

if isa(obj,'hybZono')
    [leaves] = getLeaves(obj,optSolver);
    if isempty(leaves)
        error('Empty set.')
    end
    nLeaves = size(leaves,2);
    for i = 1:nLeaves
        Zi = conZono(obj.Gc,obj.c+obj.Gb*leaves(:,i),obj.Ac,obj.b-obj.Ab*leaves(:,i));
        Ci = complement(Zi,bounds);
        if i == 1
            out = Ci;
        else
            out = out & Ci;
        end
    end
    return
end

% Determine min factor value needed for scaled obj to contain bounds interval


end