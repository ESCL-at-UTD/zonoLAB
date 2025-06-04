% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Class:
%       Abstract zonotope used as superclass for hybZono, conZono, and zono.
%   Syntax:
%       None
%   Inputs:
%       None
%   Outputs:
%       None
%   Notes:
%       Cannot be instantiated. Used to hold properties and methods shared
%       by hybZono, conZono, and zono classes.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
classdef (Abstract) abstractZono < DisplayNonScalarObjectAsTable
    
    properties (Abstract) % These properties must be defined by each subclass
        G   double  % Generator matrix (n x nG)
        Gc  double  % Continuous generator matrix (n x nGc)
        Gb  double  % Binary generator matrix (n x nGb)
        c   double  % Center (n x 1)
        A   double  % Constraint matrix (nC x nG)
        Ac  double  % Continuous constraint matrix (nC x nGc)
        Ab  double  % Binary constraint matrix (nC x nGb)
        b   double  % Constraint vector (nC x 1)
        n   int32   % Dimension
        nG  int32   % Number of generators
        nGc int32   % Number of continuous generators
        nGb int32   % Number of binary generators
        nC  int32   % Number of constraints
    end

    methods
        
        % Arithmetic
        obj = and(obj1,obj2,R)      % Generalized intersection
        obj = boundingBox(obj)      % Bounding box
        obj = cartProd(obj1,obj2)   % Cartesian product
        obj = convexHull(obj1,obj2) % Convex hull
        obj = mtimes(M,obj)         % Linear mapping
        obj = plus(obj1,obj2)       % Minkowski sum
        obj = pontryDiff(obj1,obj2)  % Pontryagin difference
        obj = projection(obj,dims)  % Projection onto specified dimensions
        [NN,Y] = reluNN(X,Ws,bs,a)  % Input-output mapping of a ReLU neural network
        out = stepMLD(X0,U,W,A,B_u,B_w,B_aff,E_x,E_u,E_w,E_aff) % 1-step reachable set for MLD
        [s,x] = supportFunc(obj,d)  % Support function
        obj = union(obj1,obj2)      % Union

        % Visualization
        [v,f] = plot(obj,varargin)  % Plot and output vertices and faces

        % Auxiliary 
        [out1, out2] = matchSetType(obj1,obj2) % Output sets as the same class (zono, conZono, or hybZono)

    end



    %% Hidden/Dependent Parameters
    properties (Hidden,Dependent)
        ub
        lb
        bounds
    end
    % Getter/Setter
    methods
        function out = get.bounds(obj)
            bb = obj.boundingBox;
            ub = bb.c + sum(abs(bb.G),2);
            lb = bb.c - sum(abs(bb.G),2);
            out = [lb,ub];
        end
        function out = get.ub(obj)
            bb = obj.boundingBox;
            out = bb.c + sum(abs(bb.G),2);
        end
        function out = get.lb(obj)
            bb = obj.boundingBox;
            out = bb.c - sum(abs(bb.G),2);
        end
    end

end