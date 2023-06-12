classdef conZono < abstractZono
    % Constrained zonotope class of the form
    %   Z = { c + G \xi | ||\xi||_inf <= 1, A \xi = b }
    % Define as:
    %   Z = conZono(G,c,A,b); G - generator matrix, c - center, 
    %                         A - constraint matrix, b - constraint vector
    %   Z = conZono(G,c); G - generator matrix, c - center, A = [], b = []
    %   Z = conZono(z); z - a zono object, A = [], b = []

    properties
        G       % Generator matrix (n x nG)
        c       % Center (n x 1)
        A       % Constraint matrix (nC x nG)
        b       % Constraint vector (nC x 1)
    end
    properties (Dependent) % These properties get automatically updates when used
        n       % Dimension
        nG      % Number of generators
        nC      % Number of constraints
    end
    properties (Hidden) % Unused properties from superclass
        Gc      % Continuous generator matrix (n x nGc)
        Gb      % Binary generator matrix (n x nGb)
        Ac      % Continuous constraint matrix (nC x nGc)
        Ab      % Binary constraint matrix (nC x nGb)
        nGc     % Number of continuous generators
        nGb     % Number of binar generators
    end

    methods
        % Constructor
        function obj = conZono(varargin)
            if nargin == 4
                obj.G = varargin{1};
                obj.c = varargin{2};
                obj.A = varargin{3};
                obj.b = varargin{4};
            elseif nargin == 2
                obj.G = varargin{1};
                obj.c = varargin{2};
                obj.A = zeros(0,obj.nG);
                obj.b = [];
            elseif nargin == 1
                switch class(varargin{1})
                    case 'zono'
                        obj.G = varargin{1}.G;
                        obj.c = varargin{1}.c;
                        obj.A = zeros(0,obj.nG);
                        obj.b = [];
                    case 'double'
                        obj.G = zeros(length(varargin{1}),0);
                        obj.c = varargin{1}.c;
                        obj.A = [];
                        obj.b = [];
                end                    
            else
                error('Incorrect number of inputs.')
            end
            try [obj.c obj.G];
            catch error('Center (c) and generator matrix (G) must have the same number of rows.')
            end
            try [obj.G; obj.A];
            catch error('Generator matrix (G) and constraint matrix (A) must have the same number of columns.')
            end
            try [obj.A obj.b];
            catch error('Constraint matrix (A) and vector (b) must have the same number of rows.')
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
        % Number of constraints
        function nConstraints = get.nC(obj)
            nConstraints = size(obj.A,1);
        end

        % Constrained zonotope-specific methods
        [v,f] = plotConZono1D(obj,optSolver);
        [v,f] = plotConZono2D(obj,optSolver);
        [v,f] = plotConZono3D(obj,optSolver);
    end
end