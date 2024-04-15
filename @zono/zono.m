% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Class:
%       Zonotope of the form:
%       Z = { c + G \xi | ||\xi||_inf <= 1 }
%   Syntax:
%       Z = zono(G,c)
%   Inputs:
%       G - n x nG matrix to define zonotope in R^n with nG generators
%       c - n x 1 vector to define center
%   Outputs:
%       Z - zonotope as zono object
%   Notes:
%       Inherits methods from abstractZono class
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
classdef zono < abstractZono
    
    properties
        G       % Generator matrix (n x nG)
        c       % Center (n x 1)
    end
    properties (Dependent) % These properties get automatically updated when used
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
        nGb     % Number of binary generators
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
            % Dimension compatibility checking
            try [obj.c obj.G];
            catch error('Center (c) and generator matrix (G) must have the same number of rows.')
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
        [v,f] = plotZono1D(obj);    % Plot in 1 dimension
        [v,f] = plotZono2D(obj);    % Plot in 2 dimensions
        [v,f] = plotZono3D(obj);    % Plot in 3 dimensions
        [out] = checkSetContain(obj,X,varargin);
        %out = nCross(H);
        %out = removeColumns(M, List);
        [C, d] = zono2HPoly(obj, varargin);
        [F_approx, s, alpha] = outerMinRPI(A, W, epsilon, max_iter, varargin)
    end
end