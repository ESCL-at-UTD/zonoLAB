% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Class:
%       User-defined solver options
%   Syntax:
%       optSolver = solverOptions('property1',value1,'property2',value2,...);
%   Inputs:
%       See property/value pairs below for available options
%   Outputs:
%       optSolver - solver options as solverOptions object
%   Notes:
%       lpSolver - only linprog and gurobi are currently supported
%       milpSolver - only gurobi is currently supported 
%       MIPFocus - 0: balance finding new feasible solutions and proving optimality
%                  1: focus on finding feasbile solutions quickly
%                  2: foucs on proving optimality
%                  3: focus on bounding objective
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
classdef solverOptions 
    % User-defined solver options.
    % Supported options:
    % Linear Programs (LPs): linprog, gurobi
    % Mixed-Integer Linear Programs (MILPs): intlinprog, gurobi
    
    properties
        lpSolver = 'gurobi'         % 'gurobi' or 'linprog'
        milpSolver = 'gurobi'       % 'gurobi' ('intlinprog' eventually)
        nSolutions = 1;             % Default number of solutions for MILP (integer between 1 and 2e9)
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