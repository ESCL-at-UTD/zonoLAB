% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return vertices and faces for a zonotopic set and plot the set 
%   Syntax:
%       [v,f] = plot(Z,optPlot)
%       [v,f] = plot(Z,color)
%       [v,f] = plot(Z,color,alpha)
%   Inputs:
%       Z - zonotopic set in 1D, 2D, or 3D (hybZono, conZono, or zono object)
%       optPlot - plotting options (plotOptions object) 
%       color - face color, specify MATLAB patch compatible color
%       alpha - face transparency, specify scalar in interval [0,1] 
%   Outputs:
%       v - nV x n matrix, each row denoting the positions of the nV vertices
%       f - nF x nMax matrix, each row denoting the vertices (up to nMax) contained
%                          in the nF faces (padded with NaN if face
%                          contains less than nMax vertices)
%   Notes:
%       Plotting of conZono and hybZono is highly dependent on tolerances
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function varargout = plot(obj,varargin)

nOutputs = nargout;
if (nOutputs == 1) || (nOutputs >= 3) 
    error('Should either have zero outputs or 2 outputs (vertices and faces).')
end
varargout = cell(1,nOutputs);

% Parse plotting options
optPlot = plotOptions;
if length(varargin) == 1
    if isa(varargin{1},'plotOptions')
        optPlot = varargin{1};
    else
        optPlot.FaceColor = varargin{1};
    end
elseif length(varargin) == 2
    optPlot.FaceColor = varargin{1};
    optPlot.FaceAlpha = varargin{2};
end
P = plotOptionsStruct(optPlot);

% Add cases for empty sets and sets with a single point %%%%%%%%%%%%%%%%%%
if (obj.n > 3)
    disp(['Can only plot in 1, 2, or 3 dimensions.'])
elseif obj.n == 1
    switch class(obj)
        case 'zono'
            [v,f] = plotZono1D(obj);
        case 'conZono'
            [v,f] = plotConZono1D(obj,optPlot.SolverOpts);
        case 'hybZono'
            [v,f] = plotHybZono1D(obj,optPlot.SolverOpts);
    end
elseif obj.n == 2
    switch class(obj)
        case 'zono'
            [v,f] = plotZono2D(obj);
        case 'conZono'
            [v,f] = plotConZono2D(obj,optPlot.SolverOpts);
        case 'hybZono'
            [v,f] = plotHybZono2D(obj,optPlot.SolverOpts);
    end
else
    switch class(obj)
        case 'zono'
            [v,f] = plotZono3D(obj);
        case 'conZono'
            [v,f] = plotConZono3D(obj,optPlot.SolverOpts);
        case 'hybZono'
            [v,f] = plotHybZono3D(obj,optPlot.SolverOpts);
    end
    view(3)
end

if strcmp(optPlot.Display,'on')
    % Plot as a single patch with a single handle
    P.Faces = f;
    P.Vertices = v;
    patch(P)
elseif strcmp(optPlot.Display,'off')
    % Do nothing
elseif strcmp(optPlot.Display,'individual')
    % Plot as multiple patches with a handle for each face
    for indx = 1:size(f,1)
        P.Faces = f(indx,:);
        P.Vertices = v;
        patch(P)
    end
end

if nOutputs == 2
    varargout{1} = v;
    varargout{2} = f;
end    

end