%	Hybrid zonotope object function: rescales the continuous generators of 
%			a hybrid zonotope to reduce error while removing constraints. 
%			This method also removes constant generators.
%			*solves 2ngc MILPs to find bounds on the factors
% 
%	hybZono/rescaleCont
% 
%	Syntax: 
%		[ Zh n, box , boxBin ] = rescaleCont(Zh)
% 
%	Inputs:
%		Zh : n dimensional hybrid zonotope in HCG-rep
% 
%	Outputs:
%		Zh : n dimensional hybrid zonotope in HCG-rep with rescaled continuous factors
% 		box : ngc x 2 dimensional vector containing bounds on each of the 
%				factors
%		boxBin : binary factors for the values obtained in box
%	
%	Trevor Bird - bird6@purdue.edu - Purdue 2022
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [ Zh , box , boxInt ] = rescaleCont(Zh,solverTol,zeroTol)
if nargin < 2
% 	solverTol = max(max(abs([Zh.Ac Zh.Ab Zh.b])))*(1e-6);
	solverTol = norm([Zh.Ac Zh.Ab Zh.b],inf)*(1e-6);
end
if nargin < 3
	zeroTol = (Zh.nGc+Zh.nGb+1)*eps('double')*norm([Zh.Ac Zh.Ab Zh.b],inf);	% this is matlabs default for rref..
end
if ~Zh.nC
	% there is nothing to do.. rescaling would be from [-1 1] -> [-1 1]
	box = [ -ones(Zh.nGc,1) , ones(Zh.nGc,1) ];
	boxInt = box;
	return
end

% define the gurobi model
% concatinate continuous and binary factors and convert from {-1 1} to {0 1} binaries
model.A = sparse([ Zh.Ac 2*Zh.Ab ]);
if isempty(Zh.Ab)
	model.rhs = Zh.b;
else
	model.rhs = Zh.b+Zh.Ab*ones(Zh.nGb,1);
end
% variable type continuous or binary
model.vtype(1:Zh.nGc) = 'C';
model.vtype(Zh.nGc+1:Zh.nGc+Zh.nGb) = 'B';
% Ax=b 
model.sense = '=';

% define the solver options to look for all feasible solutions
% params.PoolSearchMode = 0;		% 0 => find X feasible solutions
% params.MIPFocus = 1;				% 1 => just look for feasible solutions
params.OutputFlag = 0;				% 0 => don't output to console 

% ************************************* %
params.Threads = 4;	% number of cores -> 0 => will use however many it wants
% ************************************* %

% initialize the box
box = [ -1.1*ones(Zh.nGc,1) , 1.1*ones(Zh.nGc,1) ] ;

if Zh.nC
	% which generators have nonzero values in the constraint matrix?
	inds = 1:Zh.nGc;
	inds = inds(any(abs(Zh.Ac)>zeroTol));

	modelT = model;

	% only rescale factors that are present in the constraints
	for i = inds

		model = modelT;
		% search for i variable bounds
		s = zeros(Zh.nGc+Zh.nGb,1);
		s(i) = 1;
		model.obj = s;
		% add bounds for inf norm
		model.ub = ones(Zh.nGc+Zh.nGb,1);
		model.lb = [ -1*ones(Zh.nGc,1); zeros(Zh.nGb,1) ];
		% relax inf norm constraint on ith gen
		model.ub(i) = inf;
		model.lb(i) = -inf;
		
		% lower bound
		model.modelsense = 'min';
		result = gurobi(model,params);

		if isfield(result,'x')
			if isfield(result,'mipgap')
				if result.mipgap > 0
					% if the gap is not zero.. better to make this choices using upper bound
					box(i,1) = result.objbound;
% 					boxT(1) = result.objbound;
% 					box(i,1) = result.objval;
				else
					box(i,1) = result.objval;
% 					boxT(1) = result.objval;
				end
			else
				box(i,1) = result.objval;
% 				boxT(1) = result.objval;
			end
			% upper bound -> assume that if lower bound was infeasible then
			% this one is too
			model.modelsense = 'max';
			result = gurobi(model,params);
			if isfield(result,'x')
				if isfield(result,'mipgap')
					if result.mipgap > 0
						% if the gap is not zero.. better to make this choices using upper bound
						box(i,2) = result.objbound;
% 						boxT(2) = result.objbound;
	% 					box(i,2) = result.objval;
					else
						box(i,2) = result.objval;
% 						boxT(2) = result.objval;
					end
				else
					box(i,2) = result.objval;
% 					boxT(2) = result.objval;
				end
			else
				stop = 1;
	% 			error('.. why is this infeasible / unbounded?')
			end
		else
			stop = 1;
% 			error('.. why is this infeasible / unbounded?')
		end


% 		if isfield(result,'x')
% 			box(i,1) = result.x(i);
% % 			boxBin{indsCheck(i)} = 2*opt.xopt(Zh.ngc+1:end) - 1;
% 		end
		
% 		if isfield(result,'x')
% 			box(i,2) = result.x(i);
% % 			boxBin{indsCheck(i)} = 2*opt.xopt(Zh.ngc+1:end) - 1;
% 		end
% 		box{i} = boxT;
% 		box(i,:) = boxT;
	end
	
else
% 	error('There are no constraints to rescale..')
end

% saturate with inf norm constraints
boxInt = box;
boxInt(boxInt < -1) = -1;
boxInt(boxInt > 1) = 1;

% find equivalent variable in terms of these bounds
xir = (boxInt(:,2) - boxInt(:,1))/2;
xim = (boxInt(:,2) + boxInt(:,1))/2;

indsLessThanTol = xir<solverTol;
if any(indsLessThanTol)
	xir(indsLessThanTol) = 0;
% 	xim(indsLessThanTol) = 1;
end

% and rescale the continuous generators
Zh.c = Zh.c + Zh.Gc*xim;
Zh.Gc = Zh.Gc*diag(xir);
Zh.b = Zh.b - Zh.Ac*xim;
Zh.Ac = Zh.Ac*diag(xir);

% % rescale the box
% for i = 1:size(box,1)
% 	if xim(i) < 1
% 		box(i,1) = (box(i,1)-xim(i))/xir(i);
% 		box(i,2) = (box(i,2)-xim(i))/xir(i);
% 	end
% end

end