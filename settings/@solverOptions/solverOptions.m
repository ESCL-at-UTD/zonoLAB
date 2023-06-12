classdef solverOptions 
    % User-defined solver options.
    % Supported options:
    % Linear Programs (LPs): linprog, gurobi
    % Mixed-Integer Linear Programs (MILPs): intlinprog, gurobi
    
    properties
        lpSolver = 'linprog'        % Default comes with MATLAB
        milpSolver = 'intlinprog'   % Default comes with MATLAB
        nSolutions = 1;             % Default number of solutions for MILP
        MIPFocus = 0;               % Default solution strategy for MILP
    end
    
    methods
        % Constructor
        function obj = solverOptions(varargin)
            if mod(length(varargin),2) % Check if even
                error('Inputs must be specified in name-value pairs.')
            end
            for i = 1:2:(length(varargin)-1)
                obj.(varargin{i}) = varargin{i+1};
            end
        end
    end
end

