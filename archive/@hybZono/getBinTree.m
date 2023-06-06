%	Hybrid zonotope object function: branch the previous iterations binary 
%			tree with N nonempty leaves by solving N MILPs
% 
%	hybZono/getBinTree
% 
%	Syntax: 
%		Zh = getBinTree(Zh)
% 
%	Inputs:
%		Zh : n dimensional hybrid zonotope object in HCG-rep
% 
%	Outputs:
% 		Zh : n dimensional hybrid zonotope object in HCG-rep with updated 
%				binary tree
%	
% NOTE: Runs N MIPs to find feasible paths of binary tree, one for each of
% the N nonempty leaves found during previous iterations. If no binary tree
% has been previously defined, solves one large MIP to find leaves.
% 
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function Zh = getBinTree(Zh)
Zh = getDimensions(Zh);
nd = size(Zh.binTree,1);

if nd
	
	% branch from each leaf
	Tip = [];
	for i = 1:size(Zh.binTree,2)
		
		% define leaf from last iteration to search new branches for nonempty leaves
		Tim = Zh.binTree(:,i);
		branch = hybZono(Zh.Gc,Zh.Gb(:,nd+1:end),Zh.c+Zh.Gb(:,1:nd)*Tim,Zh.Ac,Zh.Ab(:,nd+1:end),Zh.b-Zh.Ab(:,1:nd)*Tim);
		
		% find the new leaves branching from this interior node
		
		Tit1 = getLeaves(branch);
		
		Tit2 = [ repmat(Tim,1,size(Tit1,2)) ; Tit1 ];
		Tip = [ Tip , Tit2 ];
		
	end

else
	
	% doesn't have a tree yet.. search whole space
	Tip = getLeaves(Zh);
	
end

Zh.binTree = Tip;

end	% end of getBinTree_MILP