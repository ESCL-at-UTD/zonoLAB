%% Plot Function for memZono
function plot(obj,dims,varargin)
    arguments
        obj memZono
        dims = [];
    end
    arguments (Repeating)
        varargin
    end
    
    if isempty(dims), dims = obj.dimKeys; end
    if length(dims) > 3, error('specify dims... too many to plot'); end

    Z_ = obj.projection(dims).Z;
    plot(Z_,varargin{:}); 
end