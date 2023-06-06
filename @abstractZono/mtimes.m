function out = mtimes(M,obj)

% Standardized header

% Check for compatable dimensions
if (~isscalar(M)) && (~ismatrix(M))
    error('Only scalar multiplication and linear transformations by matrices are supported.')
elseif (~isscalar(M)) && (size(M,2) ~= obj.n)
    error(['First argument must be a scalar or matrix with ',num2str(obj.n),' columns.'])
end

switch class(obj)
    case 'zono'
        out = zono(M*obj.G,M*obj.c);
    case 'conZono'
        out = conZono(M*obj.G,M*obj.c,obj.A,obj.b);
    case 'hybZono'
        out = hybZono(M*obj.Gc,M*obj.Gb,M*obj.c,obj.Ac,obj.Ab,obj.b);
end
end