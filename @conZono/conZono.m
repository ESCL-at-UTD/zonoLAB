% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Class:
%       Constrained zonotope of the form:
%       Z = { c + G \xi | ||\xi||_inf <= 1, A \xi = b }
%   Syntax:
%       Z = conZono(G,c,A,b)
%       Z = conZono(G,c)
%       Z = conZono(z)
%       Z = conZono([H f])
%   Inputs:
%       G - n x nG matrix to define constrained zonotope in R^n with nG generators
%       c - n x 1 vector to define center
%       A - nC x nG matrix to define nC equality constraints (A \xi = b)
%       b - nC x 1 vector to define nC equality constraints (A \xi = b)
%       z - zono object to be recast as constrained zonotope
%       H - nH x n matrix for set in H-rep (H x <= f)
%       f - nH x 1 vector for set in H-rep (H x <= f)
%   Outputs:
%       Z - constrained zonotope as conZono object
%   Notes:
%       Inherits methods from abstractZono class
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
classdef conZono < abstractZono

    properties
        G       % Generator matrix (n x nG)
        c       % Center (n x 1)
        A       % Constraint matrix (nC x nG)
        b       % Constraint vector (nC x 1)
    end
    properties (Dependent) % These properties get automatically updated when used
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
        nGb     % Number of binary generators
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
                        if size(varargin{1},2) == 1
                            obj.G = zeros(length(varargin{1}),0);
                            obj.c = varargin{1};
                            obj.A = [];
                            obj.b = [];
                        else
                            obj = hPoly2conZono(varargin{1});
                        end
                end                    
            else
                error('Incorrect number of inputs.')
            end
            % Dimension compatibility checking
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
        [v,f] = plotConZono1D(obj,optSolver);   % Plot in 1 dimension
        [v,f] = plotConZono2D(obj,optSolver);   % Plot in 2 dimensions
        [v,f] = plotConZono3D(obj,optSolver);   % Plot in 3 dimensions
        [out] = checkEmpty(obj);
        [out] = checkPointContain(obj,point);
    end
end