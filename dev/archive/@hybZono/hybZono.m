%	Hybrid zonotope object: this object contains the set definition and 
%		all code necessary for its set operations
%	
%   Constructor: 
%       Z = hybZono(Gc,Gb,c,Ac,Ab,b)
%       Outputs
%			Z  : hybrid zonotope
%		Inputs
%			-- Generators   : z = Gc*xic + Gb*xib + c --
%			Gc : (n x ngc)  continuous generator matrix 
%			Gb : (n x ngb)  binary generator matrix 
%           c  : (n x 1)    center vector 
%			-- Constraints  : Ac*xic + Ab*xib = b --
%			Ac : (nc x ngc) continuous equality constraints matrix 
%			Ab : (nc x ngb) binary equality constraints matrix 
%			b  : (nc x 1)   equality constraints vector 
%		Other Properties
%			-- Dimensions --
%			n   : dimension of space
%			ngc : number of continuous generators
%			ngb : number of binary generators
%			nc  : number of constraints
%			-- Binary tree to store feasible paths --
%			binTree : (ngb x nfeas) each column a feasible combination of binary factors
%   References:
%   Acknowledgements: 
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
classdef hybZono < DisplayNonScalarObjectAsTable
    
    properties
		% center and generators: z = Gc*xic + Gb*xib + c
		Gc = [];	% (n x ngc) matrix
		Gb = [];	% (n x ngb) matrix
        c  = [];	% (n x 1) vector
		% constraints on xi space: s.t. Ac*xic + Ab*xib = b
		Ac = [];	% (nc x ngc) matrix
		Ab = [];	% (nc x ngb) matrix
		b  = [];	% (nc x 1) vector
		% binary tree containing feasible combinations
		binTree = []; % (ngb x nfeas) matrix of 1 and -1 entries
		% dimesions of set: all scalar
		n   = [];		
		ngc = [];
		ngb = [];
		nc  = [];
	end
    
    methods
        
        %--% class constructor %--%
        function obj = hybZono(varargin)
            % construct an instance of this class
            if nargin == 1 % construct a point 
                obj.c = varargin{1};
            elseif nargin == 2 % construct a zonotope
                obj.Gc = varargin{1};
                obj.c  = varargin{2};
            elseif nargin == 3 % construct an unconstrained hybrid zonotope
                obj.Gc = varargin{1};
                obj.Gb = varargin{2};
                obj.c  = varargin{3};
            elseif nargin == 4 % construct a constrained zonotope
                obj.Gc = varargin{1};
                obj.c  = varargin{2};
                obj.Ac = varargin{3};
                obj.b  = varargin{4};
            elseif nargin == 6 % construct a hybrid zonotope
                obj.Gc = varargin{1};
                obj.Gb = varargin{2};
                obj.c  = varargin{3};
                obj.Ac = varargin{4};
                obj.Ab = varargin{5};
                obj.b  = varargin{6};
			else
               error(['Incorrect input arguments. \nAccepted inputs are: ', ...
                      '\nZ = hybZono(c), \nZ = hybZono(Gc,c), \nZ = hybZono(Gc,Gb,c), ',...
                      '\nZ = hybZono(Gc,c,Ac,b), \nZ = hybZono(Gc,Gb,c,Ac,Ab,b).'],[])
            end
            obj = getDimensions(obj);
        end
        
        %--% set operations %--%
% 		out = mtimes(R,obj)		% linear mapping
% 		out = plus(obj,W)		% Minkowski sum
% 		out = and(obj,Y,R)		% generalized intersection
%         out = cartProd(obj,Y)	% cartesian product
% 		out = zIntH(obj,h,f)	% intersection with hyperplane
% 		out = union(obj,V)		% union of two hybrid zonotopes
% 		out = pontryDiff(obj,Z2)% Pontryagin difference between hybrid zonotope and a zonotope
        
        %--% analysis and plotting %--%
        out = getDimensions(obj)			% get dimensions of the hybrid zonotope
        out = isempty(obj)					% solve MILP to find if set is empty
        out = checkZintH(obj,h,f)			% check if set Z intersects the halfspace
%         [out1,out2] = pointContain(obj,p)	% solve MILP to find if points belong to the set 
% 		out = intervalBox(obj)				% find interval box containing the set by solving 2n MILPS
        out = getLeaves(obj)				% get leaves.. will call gurobi if it exists otherwise uses LPsearch below
		out = getLeaves_gurobi(obj)			% find the integer feasible set using gurobi 
%         out = getLeaves_LPsearch(obj)		% find the integer feasible set with homemade MILP exhaustive search.. gurobi is better
        out = getBinTree(obj)				% iteratively find the integer feasible set by branching on previous leaves
        out = plotBinaryTree(obj)			% plot the binary tree for feasible combinations of binary factors
        [out1,out2] = plot_MPT(obj,color,alpha,dims)	% plot the equivalent complex of convex sets 
        out = hyb2polyhedron_MPT(obj,combos)	% convert hybrid zonotope to a complex of convex polyhedron using MPT 3.0
%         out = hyb2polyhedron_ApproxLPs(obj,res,dims,combos)	% convert hybrid zonotope to an inner approximation complex of convex polyhedron using a shooting method
    end
end

