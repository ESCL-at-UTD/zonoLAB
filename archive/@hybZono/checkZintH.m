%	Hybrid zonotope object function: check if any element in Z is in the
%									 halfspace defined by hz<=f
%	hybZono/checkZintH
% 
%	Syntax: 
%		B = checkZintH(Zh,h,f)
% 
%	Inputs:
%		Zh : n dimensional hybrid zonotope object in HCG-rep
%		h : (1 x ng) dimensional vector such that hz<=f for all z in Z
%		f : (1) scalar such that hz<=f for all z in Z
% 
%	Outputs:
% 		B : bool where 1-> hz<=f for some z in Z, 0-> otherwise
%	
%	References: 
%		[1] V. Raghuraman et al. "Set operations and order reductions for
%		constrained zonotopes" 2020
% 
% Notes: This function checks if any element in Z lies within the half
% space by solving an MILP but could be modified to detect an exact
% intersection, i.e. hz == f for some z in Z
% 
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function B = checkZintH(obj,h,f)
Zh = getDimensions(obj);

% check intersection using MPT to solve a MILP 
% variable type continuous or binary
vtype(1:Zh.ngc) = 'C';
vtype(Zh.ngc+1:(Zh.ngc+Zh.ngb)) = 'B';
% equality constraints of hybrid zono
Ae = [ Zh.Ac, 2*Zh.Ab ];
be = Zh.b+Zh.Ab*ones(Zh.ngb,1);
% inequality constraints for half space intersection
A = h*[ Zh.Gc 2*Zh.Gb ];
b = f - h*Zh.c;

% check if feasible
milp = Opt('A',A,'b',b,'Ae',Ae,'be',be,'lb',-1*ones(Zh.ngc+Zh.ngb,1),'ub',1*ones(Zh.ngc+Zh.ngb,1),'vartype',vtype); % formulates the MILP
opt = mpt_solve(milp); % solves the MILP

if opt.exitflag ~= 1
	B = false;
else
	B = true;
end

end		% end checkZintH