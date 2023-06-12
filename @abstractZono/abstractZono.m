classdef (Abstract) abstractZono < DisplayNonScalarObjectAsTable
    % This abstract class defines the properties and methods used by the
    % zono, conZono, and hybZono subclasses
    
    properties (Abstract) % These properties must be defined by each subclass
        G       % Generator matrix (n x nG)
        Gc      % Continuous generator matrix (n x nGc)
        Gb      % Binary generator matrix (n x nGb)
        c       % Center (n x 1)
        A       % Constraint matrix (nC x nG)
        Ac      % Continuous constraint matrix (nC x nGc)
        Ab      % Binary constraint matrix (nC x nGb)
        b       % Constraint vector (nC x 1)
        n       % Dimension
        nG      % Number of generators
        nGc     % Number of continuous generators
        nGb     % Number of binar generators
        nC      % Number of constraints
    end

    methods
        
        % Arithmetic
        obj = plus(obj1,obj2)       % Minkowski sum
        obj = cartProd(obj1,obj2)   % Cartesian product
        obj = mtimes(M,obj)         % Linear mapping

        % Plot
        [v,f] = plot(obj,varargin)  % Plot and output vertices and faces
    end
end