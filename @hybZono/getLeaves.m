function [leaves] = getLeaves(obj,opts)

% if there are no constraints, then the tree is full
if isempty(obj.b)
	combos = ff2n(obj.nGb);
	combos(combos==0) = -1;
	leaves = combos';
	return
end

if isempty(opts)
    opts = solverOptions;
end

if ~strcmp(opts.milpSolver,'gurobi')
    error('Gurobi is currently required for this function.')
end

opts.nSolutions = 2^obj.nGb;
opts.MIPFocus = 0;
Aeq = [obj.Ac 2*obj.Ab];
beq = obj.b+obj.Ab*ones(obj.nGb,1);
lb = -ones(obj.nGc+obj.nGb,1);
ub =  ones(obj.nGc+obj.nGb,1);
lb((obj.nGc+1):end) = 0; % Set lower bound on binary factors to zero
vType(1:obj.nGc) = 'C';
vType(obj.nGc+1:(obj.nGc+obj.nGb)) = 'B';
[x,~,~] = solveMILP([],[],[],Aeq,beq,lb,ub,vType,opts);

% convert binaries back to {-1 1}
leaves = (x((obj.nGc+1):end,:)-0.5)/0.5;

end

