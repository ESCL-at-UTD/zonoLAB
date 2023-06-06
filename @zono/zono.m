classdef zono < abstractZono
    % Zonotope class of the form
    %   Z = { c + G \xi | ||\xi||_inf <= 1 }
    % Define as:
    %   Z = zono(G,c); G - generator matrix, c - center
    
    properties
        G       % Generator matrix (n x nG)
        c       % Center (n x 1)
    end
    properties (Dependent)
        n       % Dimension
        nG      % Number of generators
    end
    properties (Hidden) % Unused properties from superclass
        Gc      % Continuous generator matrix (n x nGc)
        Gb      % Binary generator matrix (n x nGb)
        A       % Constraint matrix (nC x nG)
        Ac      % Continuous constraint matrix (nC x nGc)
        Ab      % Binary constraint matrix (nC x nGb)
        b       % Constraint vector (nC x 1)
        nGc     % Number of continuous generators
        nGb     % Number of binar generators
        nC      % Number of constraints
    end

    methods
        % Constructor
        function obj = zono(varargin)
            if nargin == 2
                obj.G = varargin{1};
                obj.c = varargin{2};
            elseif nargin == 1
                obj.G = zeros(length(varargin{1}),0);
                obj.c = varargin{1};
            else
                error('Incorrect number of inputs.')
            end
            try 
                [obj.c obj.G];
            catch
                error('Center (c) and generator matrix (G) must have the same number of rows.')
            end
            if size(obj.c,2) ~= 1
                error('Center (c) must be a n x 1 vector.')
            end
        end
        
        % Dependent property getters
        % Dimension of the space
        function dimension = get.n(obj)
           dimension = size(obj.c,1);
        end
        % Number of generators
        function nGenerators = get.nG(obj)
            nGenerators = size(obj.G,2);
        end
        
        % Zonotope-specific methods
        plot(obj,varargin);
    end
end