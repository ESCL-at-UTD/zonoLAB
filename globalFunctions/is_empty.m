%	Hybrid zonotope object function: check if the set Z is empty by solving
%		an MILP
% 
%	Syntax: 
%		B = isempty(Zh)
% 
%	Inputs:
%		Zh : n dimensional hybrid zonotope object in HCG-rep
% 
%	Outputs:
% 		B : bool where 0-> if there exists some z in Zh, 1-> otherwise
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function B = is_empty(Zh)

if Zh.nC == 0
	B = false;
	return
end

B = size(getLeaves(Zh), 1) > 0;
end