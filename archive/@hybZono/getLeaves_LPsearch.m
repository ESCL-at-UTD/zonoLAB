%	Hybrid zonotope object function: determine which combinations of binary
%			factors are feasible using homemade algorithm 
% 
%	hybZono/getLeaves_LPsearch
% 
%	Syntax: 
%		leaves = getLeaves_LPsearch(Zh)
% 
%	Inputs:
%		Zh : n dimensional hybrid zonotope object in HCG-rep
% 
%	Outputs:
% 		leaves : nb x N matrix with elements in {-1,1}. Each row is a 
%					binary factor and each column is a feasible combination
%	
% NOTE: Solves relaxation of MILP at each branch node to possibly detect an
% empty branch. If branch is is found to be nonempty, the path to the
% branch is stored and branched on in the next iteration
% 
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function leaves = getLeaves_LPsearch(Zh)

% need to check each layer for feasible paths.. solve LP to approximate if set is empty
% gets more accurate as we work our way down the tree until exact at the leaves
Tt = [];
Tp = Zh.binTree;
if isempty(Zh.binTree)
	n = 1;
	Tp(:,1) = NaN;
else
	n = size(Tp,2);
end

% branch on all the layers and check if relaxation of MILP is empty
for j = 1:Zh.ngb
	for i = 1:n
		Tit = checkChildren_LP(Zh,Tp(:,i));
		Tt = [ Tt , Tit ];
	end
	Tp = Tt;
	n = size(Tp,2);
	Tt = [];
end

leaves = Tp;

end	% end of getBinTree_LPsearch

function Ti = checkChildren_LP(Xp,Tim)
% solve 2 LPs to potentially detect an empty branch

if isnan(Tim)
	Tim = [];
end

% new potential hybrid zonotope to check for next binary is = 1
T1 = [ Tim ; 1 ];
nd = size(T1,1);

% solve LP to potentially detect empty descendents 
B1 = false;
ngc = size(Xp.Gc,2) + size(Xp.Gb,2) - nd;
nc = size(Xp.Ac,1);
s = [ zeros(1,ngc) 1 ];	% minimize the inf norm
% find the inf norm through constraints
Atop = [ eye(ngc) -ones(ngc,1) ];
Abot = [ -eye(ngc) -ones(ngc,1) ];
Ait = [ Atop ; Abot ];
bit = zeros(2*ngc,1);
% solve LP to determine if set of constraints is infeasible for this binary variable
% if so, this set is empty and the combination of binary factors is redundant
lp = Opt('f',s,'A',Ait,'b',bit,'Ae',[Xp.Ac, Xp.Ab(:,nd+1:end), zeros(nc,1)],'be',Xp.b-Xp.Ab(:,1:nd)*T1); % formulates the linear program
opt = mpt_solve(lp); % solves the LP
% check if inf norm is less than 1 or LP is infeasible
if (abs(opt.xopt(end)-1) >= 1e-10) || strcmp('infeasible',opt.how)
	B1 = true;
end

% check if binary is = -1
T2 = [ Tim ; -1 ];
B2 = false;
lp = Opt('f',s,'A',Ait,'b',bit,'Ae',[Xp.Ac, Xp.Ab(:,nd+1:end), zeros(nc,1)],'be',Xp.b-Xp.Ab(:,1:nd)*T2); % formulates the linear program
opt = mpt_solve(lp); % solves the LP
% check if inf norm is less than 1 or LP is infeasible
if (abs(opt.xopt(end)-1) >= 1e-10) || strcmp('infeasible',opt.how)
	B2 = true;
end

if ~B1
	Ti = T1;
else
	Ti = [];
end

if ~B2
	Ti = [ Ti T2 ];
end

end	% end of checkChildren