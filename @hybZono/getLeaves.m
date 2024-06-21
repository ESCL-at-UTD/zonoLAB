% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return the binary factors \xib \in {-1,1}^nGb that result in 
%       non-empty constrained zonotopes 
%   Syntax:
%       [leaves] = getLeaves(Z,optSolver)
%   Inputs:
%       Z - 1D, 2D, or 3D hybrid zonotope in HCG-Rep (hybZono object)
%       optSolver - solver options needed for mixed-integer linear propgrams
%   Outputs:
%       leaves - nGb x nLeaves matrix, each column denoting the \xi vector
%                corresponding to one of the nLeaves non-empty constrained
%                zonotopes
%   Notes:
%       If \xib denotes the i^th column of leaves, then the corresponding
%       non-empty constrained zonotope is of the form:
%       Z = { (c + Gb \xib) + Gc \xic | ||\xic||_inf <= 1, Ac \xi = b - Ab \xib }
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [leaves] = getLeaves(obj,optSolver)

% Full tree when there are no constrains (nLeaves = 2^nGb)
if isempty(obj.b)
	combos = ff2n(obj.nGb);
	combos(combos==0) = -1;
	leaves = combos';
	return
end

% Continue if set has equality constraints
if nargin < 2 || isempty(optSolver)
    optSolver = solverOptions;
end

if ~strcmp(optSolver.milpSolver,'gurobi')
    error('Gurobi is currently required for this function.')
end

% Get up to 2^nGb solutions
optSolver.nSolutions = 2^obj.nGb;
optSolver.MIPFocus = 0; % Find feasible solutions
% Problem data for mixed-integer linear program (MILP)
Aeq = [obj.Ac 2*obj.Ab];
beq = obj.b+obj.Ab*ones(obj.nGb,1);
lb = -ones(obj.nGc+obj.nGb,1);
ub =  ones(obj.nGc+obj.nGb,1);
lb((obj.nGc+1):end) = 0; % Set lower bound on binary factors to zero
vType(1:obj.nGc) = 'C';
vType(obj.nGc+1:(obj.nGc+obj.nGb)) = 'B';
[x,~,~] = solveMILP([],[],[],Aeq,beq,lb,ub,vType,optSolver);

% convert binaries back to {-1 1}
leaves = (x((obj.nGc+1):end,:)-0.5)/0.5;

end