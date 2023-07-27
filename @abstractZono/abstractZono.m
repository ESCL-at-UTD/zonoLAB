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
        nGb     % Number of binary generators
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