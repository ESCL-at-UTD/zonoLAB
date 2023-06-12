classdef plotOptions
    % User-defined plot options (subset of MATLAB's patch properties)
    
    properties
        EdgeAlpha           = 1
        EdgeColor           = [0 0 0]
        FaceAlpha           = 1
        FaceColor           = [0 0 0]
        LineStyle           = '-'
        LineWidth           = 0.5000
        Marker              = 'none'
        MarkerEdgeColor     = 'auto'
        MarkerFaceColor     = 'none'
        MarkerSize          = 6
        Display             = 'on' % Options: 'on', 'off', 'individual' 
        SolverOpts          = solverOptions;
    end
    
    methods
        % Constructor
        function obj = plotOptions(varargin)
            if mod(length(varargin),2) % Check if even
                error('Inputs must be specified in name-value pairs.')
            end
            for i = 1:2:(length(varargin)-1)
                obj.(varargin{i}) = varargin{i+1};
            end
        end
        % Store all options in a structure for calling patch
        function P = plotOptionsStruct(opts)
            propNames = properties(opts);
            for i = 1:10 %length(propNames)  Only use first 10 properties
                name = propNames(i);
                P.(name{1}) = opts.(name{1});
            end
        end
        
    end
end

