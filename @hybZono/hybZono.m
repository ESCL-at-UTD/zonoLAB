classdef hybZono < abstractZono
    % Hybrid zonotope class of the form
    %   Z = { c + Gc \xic +  Gb \xib | ||\xic||_inf <= 1 \xib \in {-1,1}^nGb, Ac \xic + Ab \xib = b }
    % Define as:
    %   Z = hybZono(Gc,Gb,c,Ac,Ab,b); Gc - continuous generator matrix, Gb - binary generator matrix, c - center, 
    %                                 Ac - continuous constraint matrix, Ab - continuous constraint matrix, b - constraint vector
    %   Z = hybZono(z); z - a zono object, Gc = G, Gb = [], Ac = Ab = [], b = []
    %   Z = hybZono(z); z - a conZono object, Gc = G, Gb = [], Ac = A, Ab = [], b = b

    properties
        Gc      % Continuous generator matrix (n x nGc)
        Gb      % Binary generator matrix (n x nGb)        
        c       % Center (n x 1)
        Ac      % Continuous constraint matrix (nC x nGc)
        Ab      % Binary constraint matrix (nC x nGb)        
        b       % Constraint vector (nC x 1)
    end
    properties (Dependent) % These properties get automatically updated when used
        n       % Dimension
        nGc     % Number of continuous generators
        nGb     % Number of binary generators
        nC      % Number of constraints
    end
    properties (Hidden) % Unused properties from superclass
        G       % Generator matrix (n x nG)
        A       % Constraint matrix (nC x nG)
        nG      % Number of generators
    end

    methods
        % Constructor
        function obj = hybZono(varargin)
            if nargin == 6
                obj.Gc = varargin{1};
                obj.Gb = varargin{2};
                obj.c  = varargin{3};
                obj.Ac = varargin{4};
                obj.Ab = varargin{5};
                obj.b  = varargin{6};
                if sum(size(obj.Gc)) == 0
                    obj.Gc = zeros(size(obj.c,1),0);
                end
                if sum(size(obj.Gb)) == 0
                    obj.Gb = zeros(size(obj.c,1),0);
                end
                if sum(size(obj.Ac)) == 0
                    obj.Ac = zeros(0,size(obj.Gc,1));
                end
                if sum(size(obj.Ab)) == 0
                    obj.Ab = zeros(0,size(obj.Gb,1));
                end
            elseif nargin == 1
                switch class(varargin{1})
                    case 'zono'
                        obj.Gc = varargin{1}.G;
                        obj.Gb = zeros(size(obj.Gc,1),0);
                        obj.c  = varargin{1}.c;
                        obj.Ac = zeros(0,size(obj.Gc,2));
                        obj.Ab = [];
                        obj.b  = [];
                    case 'conZono'
                        obj.Gc = varargin{1}.G;
                        obj.Gb = zeros(size(obj.Gc,1),0);
                        obj.c  = varargin{1}.c;
                        obj.Ac = varargin{1}.A;
                        obj.Ab = [];
                        obj.b  = varargin{1}.b;
                    case 'double'
                        obj.Gc = zeros(length(varargin{1}),0);
                        obj.Gb = zeros(length(varargin{1}),0);
                        obj.c  = varargin{1}.c;
                        obj.Ac = [];
                        obj.Ab = [];
                        obj.b  = [];
                end
            else
                error('Incorrect number of inputs.')
            end
            % Dimension compatibility checking
            try [obj.c obj.Gc obj.Gb];
            catch error('Center (c) and generator matrices (Gc and Gb) must have the same number of rows.')
            end
            try [obj.Gc; obj.Ac];
            catch error('Continuous generator matrix (Gc) and constraint matrix (Ac) must have the same number of columns.')
            end
            try [obj.Gb; obj.Ab];
            catch error('Binary generator matrix (Gb) and constraint matrix (Ab) must have the same number of columns.')
            end
            try [obj.Ac obj.Ab obj.b];
            catch error('Constraint matrices (Ac and Ab) and vector (b) must have the same number of rows.')
            end
        end
        
        % Dependent property getters
        % Dimension of the space
        function dimension = get.n(obj)
            dimension = size(obj.c,1);
        end
        % Number of continuous generators
        function nContinuousGenerators = get.nGc(obj)
            nContinuousGenerators = size(obj.Gc,2);
        end
        % Number of binary generators
        function nBinaryGenerators = get.nGb(obj)
            nBinaryGenerators = size(obj.Gb,2);
        end
        % Number of constraints
        function nConstraints = get.nC(obj)
            nConstraints = size(obj.Ac,1);
        end

        % Hybrid zonotope-specific methods
        [v,f] = plotAsConZono(obj,optSolver);   % Plot as the union of constrained zonotopes
        [v,f] = plotHybZono1D(obj,optSolver);   % Plot in 1 dimension
        [v,f] = plotHybZono2D(obj,optSolver);   % Plot in 2 dimensions
        [v,f] = plotHybZono3D(obj,optSolver);   % Plot in 3 dimensions
        [leaves] = getLeaves(obj,optSolver)     % Identify non-empty combinations of binary factors
        
    end
end