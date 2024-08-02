%	Hybrid zonotope object function: find and remove redundant generators /
%		constraints based on implicit bounds on the factors 
% 
%		** solves 2*ngc MILPs to determine bounds on continuous generators
% 
%	hybZono/removeRedundancy6
% 
%	Syntax: 
%		Zh = removeRedundancy6(Zh)
% 
%	Inputs:
%		Zh : n dimensional hybrid zonotope in HCG-rep
% 
%	Outputs:
% 		Zh : n dimensional hybrid zonotope in HCG-rep with redundant
%					equality constraints removed
%	
%	Trevor Bird - bird6@purdue.edu - Purdue 2022
%
%   NOT IN ZONOLAB (YET)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function Zh = removeRedundancy7(Zh)
%Zh = getDimensions(Zh);

% if there are no constraints, we don't have anything to do
if Zh.nC == 0
	return
end

% tolerances on what we consider 0
% zeroTol = 1000*(Zh.ng+1)*eps('double')*norm([Zh.Ac Zh.Ab Zh.b],inf);	% this is bigger than matlabs default for rref..
zeroTol = (Zh.nGc+Zh.nGb+1)*eps('double')*norm([Zh.Ac Zh.Ab Zh.b],inf);	% this is matlabs default for rref..
solverTol = norm([Zh.Ac Zh.Ab Zh.b],inf)*(1e-6);

% remove zeros
Zh = removeZeros(Zh,zeroTol);

stop = is_empty(Zh);
if stop
	stop = 1;
end

%% %--% precondition the constraints %--% 
% This seems to mess things up for some reason.. better to use QR to find dependent equations

% % use rref with full pivoting and normalized by rows inf norm
% Zh = rrefCondition(Zh,zeroTol);

% set solver tolerance to optimality tolerance of gurobi multiplied by inf norm of the constraints
% solverTol = norm([Zh.Ac Zh.Ab Zh.b],inf)*(1e-6);
% % solverTol = max(abs([Zh.Ac Zh.Ab Zh.b]),[],'all')*(1e-6);

%% %--% update the binary tree and remove redundancy %--%

% update the binary tree if necessary
solverOpts = solverOptions();
binTree = getLeaves(Zh, solverOpts);

% is there any redundancy in the binary factors?
if size(binTree,2) > 1
	% use QR to find linear independence
	[~, R, E] = qr(binTree',0); 
	if ~isvector(R)
		diagr = abs(diag(R));
	else
		diagr = abs(R(1));   
	end
	%Rank estimation
	r = find(diagr >= zeroTol*diagr(1), 1, 'last'); %rank estimation
	idx = E(1:r);	% I don't sorting matters..
% 	idx = sort(E(1:r));
	linIndep = binTree(idx,:);
	RT = r;
else
	RT = 1;
	linIndep = 1;
end

if RT <= Zh.nGb
	
	% now transform to lin independent binary factors
	M = binTree*pinv(linIndep);
	binTree = linIndep;
	if size(M,1) == size(Zh.Ab,2)
		Zh.Ab = Zh.Ab*M;
		Zh.Gb = Zh.Gb*M;
	else
		warning('Something went wrong with pseudo inverse..')
	end
	
	% is one of the remaining binary factors constant?
	constantBin = binTree(:,1) == binTree;
	constantBin = sum(constantBin,2);
	constantBin = constantBin == size(binTree,2);
	% subtract constant off c and b
	if any(constantBin)
		ind = find(constantBin == true);
		% find shifts in c and b from constant binary values
		Gbshift = Zh.Gb(:,ind)*binTree(ind,1);
		Abshift = Zh.Ab(:,ind)*binTree(ind,1);
		Zh.c = Zh.c + Gbshift;
		Zh.b = Zh.b - Abshift;
		% and remove the constant factor
		Zh.Ab(:,ind) = [];
		Zh.Gb(:,ind) = [];
		binTree(ind,:) = [];
	end

end


stop = is_empty(Zh);
if stop
	stop = 1;
end

%% %--% rescale the continuous generators %--%

[ Zh , box , boxInt ] = rescaleCont(Zh,solverTol);


stop = is_empty(Zh);
if stop
	stop = 1;
end

%--% if there were constants.. remove the zeros resulting from rescaling  %--%
rowZero = find(all(abs([ Zh.Ac Zh.Ab Zh.b ])<=zeroTol,2));
colZeroC = find(all(abs([ Zh.Gc ; Zh.Ac  ])<=zeroTol,1));
colZeroB = find(all(abs([ Zh.Gb ; Zh.Ab  ])<=zeroTol,1));

% remove empty constraints
Zh.Ac(rowZero,:) = [];
if ~isempty(Zh.Ab)
	Zh.Ab(rowZero,:) = [];
end
Zh.b(rowZero) = [];

% remove zero columns
Zh.Gc(:,colZeroC) = [];
if ~isempty(Zh.Ac)
	Zh.Ac(:,colZeroC) = [];
end
Zh.Gb(:,colZeroB) = [];
if ~isempty(Zh.Ab)
	Zh.Ab(:,colZeroB) = [];
end

% remove zeros from box
box(colZeroC,:) = [];
boxInt(colZeroC,:) = [];

% is there a redundant constraint identified?
% which of the inf norm constraints are redundant?
LHS = abs(abs(box(:,1))) <= 1 + solverTol;
RHS = abs(abs(box(:,2))) <= 1 + solverTol;
inds = 1:Zh.nGc;
indsRed = inds(LHS&RHS);


stop = is_empty(Zh);
if stop
	stop = 1;
end

%% check for linearly dependent rows in the equality constraints
if Zh.nC
% use QR instead of rref
Atemp = [ Zh.Ac Zh.Ab Zh.b ];

% % normalize the rows
% maxAtemp = max(abs(Atemp),[],2);
% Atemp = Atemp./maxAtemp;
% Zh.Ac = Zh.Ac./maxAtemp;
% Zh.Ab = Zh.Ab./maxAtemp;
% Zh.b = Zh.b./maxAtemp;

[~, R, E] = qr(Atemp',0); 

if ~isvector(R)
	diagr = abs(diag(R));
else
	diagr = abs(R(1));   
end
% rank estimation
r = find(diagr >= zeroTol*diagr(1), 1, 'last'); %rank estimation
% idx = sort(E(1:r));
% dont need to sort here I dont think..
idx = E(1:r);

if r < Zh.nC
	Zh.Ac = Zh.Ac(idx,:);
	Zh.Ab = Zh.Ab(idx,:);
	Zh.b = Zh.b(idx,:);
	% remove zeros
	ncOld = numel(idx);
	ngOld = size(Zh.Ac,2) + size(Zh.Ab,2);

	% remove zeros
	Zh = removeZeros(Zh,zeroTol);
	
	ncNew = Zh.nC;
	ngNew = Zh.nGb+Zh.nGc;
	
	if ncOld > ncNew || ngOld > ngNew
		stop = 1;
	end

end
end
stop = is_empty(Zh);
if stop
	stop = 1;
end

%% remove redundant equality constraints based on activation

% this totally happens in the MPC problems some times

% feasTol = 1e-6;
% suppA = sum(abs([Zh.Ac Zh.Ab ]),2);
% checkFeas1 = suppA >= Zh.b - feasTol;
% checkFeas2 = suppA <= Zh.b + feasTol;
% 
% if any(checkFeas1+checkFeas2==2)
% 	error('will this ever happen?')
% end


%% remove redundant factor identified from previous step

if ~isempty(indsRed)
	
	j = indsRed(1);

% 	% find coefficient in constraints not equal to 0
% 	rowInd = find(abs(Zh.Ac(:,j))>zeroTol,1);

% 	% find coefficient for this factor in constraint that is nearest mean
% 	rowInd = find(abs(Zh.Ac(:,j))>zeroTol);
% 	[~,rowIndMean] = min(abs(abs(Zh.Ac(rowInd,j))-Amean));
% 	rowInd = rowInd(rowIndMean);
	
	% max should be around 1 because of
	% conditioning
	[~,rowInd] = max(abs(Zh.Ac(:,j)));	
	Ajt = Zh.Ac(rowInd,j);
	
	% what if this value is too small?
	if isempty(rowInd) || abs(Zh.Ac(rowInd,j)) <= zeroTol
		error('Why did we get to this point of zeros?')
		for j = indsRed(2:end)
			rowInd = find(abs(Zh.Ac(:,j))>zeroTol);
			[~,rowIndMean] = min(abs(abs(Zh.Ac(rowInd,j))-Amean));
			rowInd = rowInd(rowIndMean);
			if ~isempty(rowInd) && abs(Zh.Ac(rowInd,j)) > zeroTol
				break
			end
		end
	end
	
	% solve equality constraints for this variable
	Ej1c = zeros(Zh.nGc,Zh.nC);
	Ej1c(j,rowInd) = 1;
	delG = Zh.Gc*Ej1c/Zh.Ac(rowInd,j);
	delA = Zh.Ac*Ej1c/Zh.Ac(rowInd,j);
	
	% and substitute it in other equations
	Gc = Zh.Gc - delG*Zh.Ac;
	c = Zh.c + delG*Zh.b;
	Ac = Zh.Ac - delA*Zh.Ac;
	b = Zh.b - delA*Zh.b;
	if ~isempty(Zh.Gb)
		Gb = Zh.Gb - delG*Zh.Ab;
		Ab = Zh.Ab - delA*Zh.Ab;
		Ab(rowInd,:) = [];
	else
		Gb = [];
		Ab = [];
	end

	% remove zero rows and columns
	Gc(:,j) = [];
	Ac(:,j) = [];
	Ac(rowInd,:) = [];
	b(rowInd) = [];

	% update set
	Zh.c = c;
	Zh.Gc = Gc;
	Zh.Gb = Gb;
	Zh.Ac = Ac;
	Zh.Ab = Ab;
	Zh.b = b;

	% if we don't have binary generators anymore reset binTree
	if ~Zh.nGb
		binTree = [];
	end

	if rank(binTree) < Zh.nGb
% 		error('This shouldnt happen here..')
	end
	
else
	% there weren't any redundant generators during rescaling.. 
	% no need to loop through them invidually
	return
end

% was the only redundant one removed already?
if numel(indsRed) == 1
	return
end

stop = is_empty(Zh);
if stop
	stop = 1;
end

%% loop through and find bounds to determine remaining redundant inf norm constraints

% initialize the box -> s.t. xi_j \in box_j for jth inf norm constraint relaxed
box = [ -1.1*ones(Zh.nGc,1) , 1.1*ones(Zh.nGc,1) ];		% initialize a bit bigger as default, i.e. not redundant

% we are going to check all remaining generators for redundancy 
j = 1;	% start at the first one..

while j <= Zh.nGc

	% is the max in this row greater than our zero tolerance?
	maxj = max(abs(Zh.Ac(:,j)));
	if maxj <= zeroTol
		% if not.. there is no reason to find its bounds
		j = j + 1;
	else
		
		% define the gurobi model
		% concatinate continuous and binary factors and convert from {-1 1} to {0 1} binaries
		model.A = sparse([ Zh.Ac 2*Zh.Ab ]);
		model.rhs = Zh.b+Zh.Ab*ones(Zh.nGb,1);
		% variable type continuous or binary
		model.vtype = 'C';
		model.vtype(1:Zh.nGc) = 'C';
		model.vtype(Zh.nGc+1:Zh.nGc+Zh.nGb) = 'B';
		% Ax=b 
		model.sense = '=';

		% define the solver options to look for all feasible solutions
		% params.PoolSearchMode = 0;		% 0 => find X feasible solutions
		params.OutputFlag = 0;				% 0 => don't output to console 
		
		% ************************************* %
		params.Threads = 4;	% number of cores -> 0 => will use however many it wants
		% ************************************* %
		
		% search for jth variable bounds
		s = zeros(Zh.nGc+Zh.nGb,1);
		s(j) = 1;
		model.obj = s;
		% add bounds for inf norm
		model.ub = ones(Zh.nGc+Zh.nGb,1);
		model.lb = [ -1*ones(Zh.nGc,1); zeros(Zh.nGb,1) ];
		% relax inf norm constraint on jth gen
		model.ub(j) = inf;
		model.lb(j) = -inf;
		
		% lower bound
		model.modelsense = 'min';
		result = gurobi(model,params);
		if isfield(result,'x')
			if isfield(result,'mipgap')
				box(j,1) = result.objbound;
			else
				box(j,1) = result.objval;
			end

			% if it was feasible, then it is worth checking the upper bound
			% upper bound
			model.modelsense = 'max';
			result = gurobi(model,params);
			if isfield(result,'x')
				if isfield(result,'mipgap')
					box(j,2) = result.objbound;
				else
					box(j,2) = result.objval;
				end
			else
	% 			error('Why is this infeasible / unbounded?')
				stop = 1;
			end
		else
% 			error('.. why is this infeasible / unbounded?')
			stop = 1;
		end

		if (box(j,1) >= -(1+solverTol)) && (box(j,2) <= (1+solverTol))
			
% 			% find coefficient in constraints not equal to 0
% 			rowInd = find(abs(Zh.Ac(:,j))>zeroTol,1);

% 			% find coefficient for this factor in constraint that is nearest mean
% 			rowInd = find(abs(Zh.Ac(:,j))>zeroTol);
% 			[~,rowIndMean] = min(abs(abs(Zh.Ac(rowInd,j))-Amean));
% 			rowInd = rowInd(rowIndMean);
			
			% should be close to 1 because of preconditioning
			[~,rowInd] = max(abs(Zh.Ac(:,j)));	
			Ajt = Zh.Ac(rowInd,j);

			% what if this value is too small?
			if isempty(rowInd) || abs(Zh.Ac(rowInd,j)) <= zeroTol
				error('Why did we get to this point of zeros?')
% 				for j = indsRed(2:end)
% 					rowInd = find(abs(Zh.Ac(:,j))>zeroTol);
% 					[~,rowIndMean] = min(abs(abs(Zh.Ac(rowInd,j))-Amean));
% 					rowInd = rowInd(rowIndMean);
% 					if ~isempty(rowInd) && abs(Zh.Ac(rowInd,j)) > zeroTol
% 						break
% 					end
% 				end
			end
			
			% solve for redundant variable and substitute it out
			Ej1c = zeros(Zh.nGc,Zh.nC);
			Ej1c(j,rowInd) = 1;
			delG = Zh.Gc*Ej1c/Zh.Ac(rowInd,j);
			delA = Zh.Ac*Ej1c/Zh.Ac(rowInd,j);

			Gc = Zh.Gc - delG*Zh.Ac;
			c = Zh.c + delG*Zh.b;

			Ac = Zh.Ac - delA*Zh.Ac;
			b = Zh.b - delA*Zh.b;
			
			if ~isempty(Zh.Gb)
				Gb = Zh.Gb - delG*Zh.Ab;
				Ab = Zh.Ab - delA*Zh.Ab;
				Ab(rowInd,:) = [];
			else
				Gb = [];
				Ab = [];
			end
			
			% remove zero rows and columns
			Gc(:,j) = [];
			Ac(:,j) = [];
			Ac(rowInd,:) = [];
			b(rowInd) = [];
			
			% update set
			Zh.c = c;
			Zh.Gc = Gc;
			Zh.Gb = Gb;
			Zh.Ac = Ac;
			Zh.Ab = Ab;
			Zh.b = b;

			stop = is_empty(Zh);
			if stop
				stop = 1;
			end

			% use rref with full pivoting and normalized by rows inf norm
% 			Zh = rrefCondition(Zh,zeroTol);
			
			% remove entry from box
			box(j,:) = [];
		else
			% move on to the next variable
			j = j + 1;
		end
	end
end

% Zh = rrefCondition(Zh,zeroTol);
% remove any left over zeros..
Zh = removeZeros(Zh,zeroTol);

stop = is_empty(Zh);
if stop
	stop = 1;
end

%% check for linearly dependent rows in the eequality constraints

% if Zh.nc
% 	% use QR instead of rref
% 	Atemp = [ Zh.Ac Zh.Ab Zh.b ];
% 	[~, R, E] = qr(Atemp',0); 
% 	
% 	if ~isvector(R)
% 		diagr = abs(diag(R));
% 	else
% 		diagr = abs(R(1));   
% 	end
% 	% rank estimation
% 	r = find(diagr >= zeroTol*diagr(1), 1, 'last'); %rank estimation
% 	idx = sort(E(1:r));
% 	RT = r;
% 	
% 	if RT < Zh.nc
% 		Zh.Ac = Zh.Ac(idx,:);
% 		Zh.Ab = Zh.Ab(idx,:);
% 		Zh.b = Zh.b(idx,:);
% 		% remove zeros
% 		Zh = removeZeros(Zh,zeroTol);
% 	end
% end

end
