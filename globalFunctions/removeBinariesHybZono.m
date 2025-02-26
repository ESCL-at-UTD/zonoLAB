function Zout = removeBinariesHybZono(Z, varargin)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
    
    if length(varargin)==2 && isequal(size(varargin{1}), [1 1]) && isa(varargin{1}, 'double')
        switch varargin{2}
            case 'first'
                idxs = 1:varargin{1};
            case 'last'
                idxs = Z.nGb - varargin{1} + 1 : Z.nGb;
            otherwise
                error('Incorrect input options to removeBinariesHybZono')
        end
    elseif length(varargin)==1 && isa(varargin{1}, 'double')
        idxs = varargin{1};
    else
        error('Incorrect input options to removeBinariesHybZono')
    end

    Zout = sharpHybZono(Z, idxs);
    Zout.Gc = [Zout.Gc Zout.Gb(:,idxs)];
    Zout.Gb(:,idxs) = [];
    Zout.Ac = [Zout.Ac Zout.Ab(:,idxs)];
    Zout.Ab(:,idxs) = [];
end