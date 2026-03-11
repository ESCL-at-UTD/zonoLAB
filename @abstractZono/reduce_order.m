
function obj = reduce_order(obj,technique,varargin)

if nargin < 2 || isempty(technique)
    technique = 'Exact';
end

switch technique    
    case 'Exact'
        obj = reduceExact(obj);

    case 'Outer'
        obj = reduceOuter(obj,varargin{:});

    case 'Inner'
        obj = reduceInner(obj,varargin{:});

    otherwise
        error('Unknown reduction technique')

end
end