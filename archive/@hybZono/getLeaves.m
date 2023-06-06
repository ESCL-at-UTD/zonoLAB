%	Hybrid zonotope object function: Solve a MIP to find feasible paths of
%			hybrid zonotopes binary tree
% 
%	hybZono/getLeaves
% 
%	Syntax: 
%		leaves = getLeaves(Z)
% 
%	Inputs:
%		Z : n dimensional hybrid zonotope object in HCG-rep
% 
%	Outputs:
% 		leaves : feasible combinations of the binary factors --> paths 
%					through binary tree to nonempty leaves
%	
% NOTE: Will try to use gurobi to solve the problem. If Gurobi is not
% installed, this function will throw warning and use homemade search
% algorithm. Can alter this code to be less annoying..
% 
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function leaves = getLeaves(Zh)
Zh = getDimensions(Zh);

% try
	leaves = getLeaves_gurobi(Zh);
% catch
% 	warning('Using homemade MILP exhaustive search.. installing gurobi will make this code run much faster. Or can write own function with equivalent MILP solver that will search entire integer feasible space.')
% 	leaves = getLeaves_LPsearch(Zh);
% end

end	% end of getLeaves