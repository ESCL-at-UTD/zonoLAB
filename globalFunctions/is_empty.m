%	Hybrid zonotope object function: check if the set Z is empty by solving
%		an MILP
% 
%	hybZono/isempty
% 
%	Syntax: 
%		B = isempty(Z)
% 
%	Inputs:
%		Zh : n dimensional hybrid zonotope object in HCG-rep
% 
%	Outputs:
% 		B : bool where 0-> if there exists some z in Zh, 1-> otherwise
% 
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function B = is_empty(Zh)

if Zh.nC == 0
	B = false;
	return
end

% define the gurobi model
% concatinate continuous and binary factors and convert from {-1 1} to {0 1} binaries
model.A = sparse([ Zh.Ac 2*Zh.Ab ]);
model.rhs = Zh.b+Zh.Ab*ones(Zh.nGb,1);
% add bounds for inf norm
model.ub = ones(Zh.nGc+Zh.nGb,1);
% model.lb = -1*model.ub;
model.lb = [ -1*ones(Zh.nGc,1) ; zeros(Zh.nGb,1) ];
% variable type continuous or binary
model.vtype(1:Zh.nGc) = 'C';
model.vtype(Zh.nGc+1:Zh.nGc+Zh.nGb) = 'B';
% Ax=b 
model.sense = '=';

% define the solver options to look for all feasible solutions
% params.PoolSearchMode =0;			% 0 => find X feasible solutions
params.MIPFocus = 1;				% 1 => just look for feasible solutions
params.OutputFlag = 0;				% 0 => don't output to console 

% ************************************* %
params.Threads = 4;	% number of cores -> 0 => will use however many it wants
% ************************************* %

% %--% parameters from gurobi tuning tool %--%
% params.SimplexPricing = 2;
% params.BranchDir = -1;
% params.DegenMoves = 0;
% params.Heuristics = 0;
% params.GomoryPasses = 0;
% params.ZeroHalfCuts = 0;
% params.PrePasses = 1;
% 
% %--% parameters from gurobi tuning tool for MIPFocus = 1 and Poolsearch = 2 %--%
% params.Method = 0;
% params.DegenMoves = 0;
% params.Heuristics = 0;
% params.StartNodeLimit = 0;
% params.VarBranch = 1;
% params.Cuts = 0;
% params.AggFill = 10;
% params.PrePasses = 1;

% %--% parameters from gurobi tuning tool for MIPFocus = 1 and Poolsearch = 2 %--%
% % ran for 3000 sec
% params.Method = 0;
% params.DegenMoves = 0;
% params.Heuristics = 0;
% params.VarBranch = 1;
% params.Cuts = 0;
% params.AggFill = 0;
% params.PrePasses = 1;

% gurobi_write(model,'HeatedRooms_12.mps',params)

% search the binary tree of the MIP
result = gurobi(model,params);

if isfield(result,'x')
	B = false;
else
	B = true;
end

end		% end isempty

% using MPT instead..
% % check intersection using MPT to solve a MILP 
% % variable type continuous or binary
% vtype(1:Zh.ngc) = 'C';
% vtype(Zh.ngc+1:Zh.ng) = 'B';
% % equality constraints of hybrid zono
% Ae = [ Zh.Ac, 2*Zh.Ab ];
% be = Zh.b+Zh.Ab*ones(Zh.ngb,1);
% 
% % check if feasible
% milp = Opt('Ae',Ae,'be',be,'lb',-1*ones(Zh.ng,1),'ub',1*ones(Zh.ng,1),'vartype',vtype); % formulates the MILP
% opt = mpt_solve(milp); % solves the MILP
% 
% if opt.exitflag ~= 1
% 	B = true;
% else
% 	B = false;
% end