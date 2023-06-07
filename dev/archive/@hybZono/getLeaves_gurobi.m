%	Hybrid zonotope object function: Solve a MIP to find feasible paths of
%			hybrid zonotopes binary tree
% 
%	hybZono/getLeaves_gurobi
% 
%	Syntax: 
%		leaves = getLeaves_gurobi(Z)
% 
%	Inputs:
%		Z : n dimensional hybrid zonotope object in HCG-rep
% 
%	Outputs:
% 		leaves : feasible combinations of the binary factors --> paths 
%					through binary tree to nonempty leaves
%	
% NOTE: Requires GUROBI. Cannot be reformulated to only require YALMIP or
% MPT because both of these algorithms have limited options and add
% variables that I do not know how to decifer. 
% 
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function leaves = getLeaves_gurobi(Zh)

% if there are no constraints, then the tree is full
if isempty(Zh.b)
	combos = dec2bin(2^Zh.ngb-1:-1:0)-'0';
	combos(combos==0) = -1;
	leaves = combos';
	return
end

% define the gurobi model
% concatinate continuous and binary factors and convert from {-1 1} to {0 1} binaries
model.A = sparse([ Zh.Ac 2*Zh.Ab ]);
model.rhs = Zh.b+Zh.Ab*ones(Zh.ngb,1);
% add bounds for inf norm
model.ub = ones(Zh.ngc+Zh.ngb,1);
model.lb = -1*model.ub;
% variable type continuous or binary
model.vtype(1:Zh.ngc) = 'C';
model.vtype(Zh.ngc+1:(Zh.ngc+Zh.ngb)) = 'B';
% Ax=b 
model.sense = '=';

% define the solver options to look for all feasible solutions
params.PoolSearchMode = 2;			% find X many feasible solutions
if 2^Zh.ngb > 2000000000
	params.PoolSolutions = 2000000000;	% limit on how many gurobi will pool
else
	params.PoolSolutions = 2^Zh.ngb;	% find a max of 2^nb solutions
end
params.MIPFocus = 1;				% just look for feasible solutions
params.OutputFlag = 0;				% don't output to console 

% ************************************* %
params.Threads = 0;	% number of cores -> default 0 will use however many it wants
% ************************************* %

% search the binary tree of the MIP
result = gurobi(model,params);

if ~isfield(result,'pool')
	warning('This hybrid zonotope is empty.. MIP detected 0 feasible solutions')
	leaves = NaN;
	return
end

% and extract the results
nl = numel(result.pool);
leaves = zeros(Zh.ngb,nl);
for i = 1:nl
	leaves(:,i) = result.pool(i).xn(Zh.ngc+1:end);
end
% convert binaries back to {-1 1}
leaves = (leaves-0.5)/0.5;

end	% end of getLeaves