%% Plot Function for memZono
% - This callse the appropriate plot method for the 
%   base class with a specification of specific dimensions
function plot(obj,dims,varargin)
    arguments
        obj memZono
        dims = [];
    end
    arguments (Repeating)
        varargin
    end
    
    if isempty(dims), dims = obj.dimKeys;
    elseif isnumeric(dims), dims = obj.dimKeys(dims);
    end

    Z_ = obj.projection(dims).Z;
    if Z_.n > 3, error('specify dims... too many to plot'); end
    plot(Z_,varargin{:});
end