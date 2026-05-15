function obj = reduceExact(obj)
% REDUCEEXACT Removes zero generators from a zonotope
%
% Input:
%   obj - zonotope object
%
% Output:
%   obj - zonotope with zero generators removed

    G = obj.G;
    c = obj.c;

    % keep only generators that are not zero
    mask = any(G ~= 0,1);
    G_reduced = G(:,mask);

    % construct reduced zonotope
    obj = zono(G_reduced,c);

end