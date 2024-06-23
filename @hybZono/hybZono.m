% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Class:
%       Hybrid zonotope of the form:
%       Z = { c + Gc \xic +  Gb \xib | ||\xic||_inf <= 1, \xib \in {-1,1}^nGb, Ac \xic + Ab \xib = b }
%   Syntax:
%       Z = hybZono(Gc,Gb,c,Ac,Ab,b)
%       Z = hybZono(z)
%       Z = hybZono(H_collection)
%       Z = hybZono({V,M})
%   Inputs:
%       Gc - n x nGc matrix to define hybrid zonotope in R^n with nGc continuous generators
%       Gb - n x nGb matrix to define hybrid zonotope in R^n with nGb binary generators
%       c  - n x 1 vector to define center
%       Ac - nC x nGc matrix to define nC equality constraints (Ac \xic + Ab \xib = b)
%       Ab - nC x nGb matrix to define nC equality constraints (Ac \xic + Ab \xib = b)
%       b  - nC x 1 vector to define nC equality constraints (Ac \xic + Ab \xib = b)
%       z - zono or conZono object to be recast as hybrid zonotope
%       H_collection - N x 1 cell array with i^th cell [H_i f_i] defining H-rep set H_i x <= f_i
%       V - n x nV matrix with each column defining a vertex in R^n
%       M - nV x N matrix indicating which of the N sets the nV vertices belong
%   Outputs:
%       Z - hybrid zonotope as hybZono object
%   Notes:
%       Inherits methods from abstractZono class
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
classdef hybZono < abstractZono

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
                    obj.Ac = zeros(0,size(obj.Gc,2));
                end
                if sum(size(obj.Ab)) == 0
                    obj.Ab = zeros(0,size(obj.Gb,2));
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
                        obj.Ab = zeros(size(obj.Ac,1),0);
                        obj.b  = varargin{1}.b;
                    case 'double'
                        obj.Gc = zeros(length(varargin{1}),0);
                        obj.Gb = zeros(length(varargin{1}),0);
                        obj.c  = varargin{1};
                        obj.Ac = [];
                        obj.Ab = [];
                        obj.b  = [];
                    case 'cell'
                        if size(varargin{1},2) == 1
                            obj = hPoly2hybZono(varargin{1});
                        else
                            obj = vPoly2hybZono(varargin{1});
                        end
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