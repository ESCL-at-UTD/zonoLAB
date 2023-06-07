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

function B = isempty(Zh)
Zh = getDimensions(Zh);

% check intersection using MPT to solve a MILP 
% variable type continuous or binary
vtype(1:Zh.ngc) = 'C';
vtype(Zh.ngc+1:(Zh.ngc+Zh.ngb)) = 'B';
% equality constraints of hybrid zono
Ae = [ Zh.Ac, 2*Zh.Ab ];
be = Zh.b+Zh.Ab*ones(Zh.ngb,1);

% check if feasible
milp = Opt('Ae',Ae,'be',be,'lb',-1*ones(Zh.ngc+Zh.ngb,1),'ub',1*ones(Zh.ngc+Zh.ngb,1),'vartype',vtype); % formulates the MILP
opt = mpt_solve(milp); % solves the MILP

if opt.exitflag ~= 1
	B = true;
else
	B = false;
end

end		% end isempty