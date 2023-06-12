function [v,f] = plot(obj,varargin)

% Standardized header

opts = plotOptions;
if length(varargin) == 1
    if isa(varargin{1},'plotOptions')
        opts = varargin{1};
    else
        opts.FaceColor = varargin{1};
    end
elseif length(varargin) == 2
    opts.FaceColor = varargin{1};
    opts.FaceAlpha = varargin{2};
end
P = plotOptionsStruct(opts);
% Add cases for empty sets and sets with a single point %%%%%%%%%%%%%%%%%%
if (obj.n > 3)
    disp(['Can only plot in 1, 2, or 3 dimensions.'])
elseif obj.n == 1
    switch class(obj)
        case 'zono'
            [v,f] = plotZono1D(obj);
        case 'conZono'
            [v,f] = plotConZono1D(obj,opts.SolverOpts);
        case 'hybZono'
            [v,f] = plotHybZono1D(obj,opts.SolverOpts);
    end
elseif obj.n == 2
    switch class(obj)
        case 'zono'
            [v,f] = plotZono2D(obj);
        case 'conZono'
            [v,f] = plotConZono2D(obj,opts.SolverOpts);
        case 'hybZono'
            [v,f] = plotHybZono2D(obj,opts.SolverOpts);
    end
else
    switch class(obj)
        case 'zono'
            [v,f] = plotZono3D(obj);
        case 'conZono'
            [v,f] = plotConZono3D(obj,opts.SolverOpts);
        case 'hybZono'
            [v,f] = plotHybZono3D(obj,opts.SolverOpts);
    end
    view(3)
end

if strcmp(opts.Display,'on')
    % Plot as a single patch with a single handle
    P.Faces = f;
    P.Vertices = v;
    patch(P)
elseif strcmp(opts.Display,'off')
    % Do nothing
elseif strcmp(opts.Display,'individual')
    % Plot as multiple patches with a handle for each face
    for indx = 1:size(f,1)
        P.Faces = f(indx,:);
        P.Vertices = v;
        patch(P)
    end
end

