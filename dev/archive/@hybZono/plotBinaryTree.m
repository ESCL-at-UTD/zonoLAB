%	Hybrid zonotope object function: plot the binary tree of Zh
% 
%	hybZono/plotBinaryTree
% 
%	Syntax: 
%		p = plotBinaryTree(Z)
% 
%	Inputs:
%		Z : n dimensional hybrid zonotope object in HCG-rep
% 
%	Outputs:
% 		p : plot handles
% 
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function p = plotBinaryTree(obj)
Xp = obj;
figure; hold on;

% first find the feasible combination of binaries
combos = getLeaves(Xp);
sumC = 0;
maxSumC = sumC;
minSumC = sumC;
p(1) = plot(0,0,'c.','markersize',30);
for j = 1:size(combos,2)
	sumC = 0;
	root = [ 0 0 ];
	incIng = size(combos,1):-1:1;
	for k = 1:size(combos,1)
		sumC = sumC + combos(k,j)*2^(incIng(k)-1);
		maxSumC = max(maxSumC,sumC);
		minSumC = min(minSumC,sumC);
		node = [ sumC , k ];
		if combos(k,j) == 1
			p(4) = plot([root(1) node(1)],[root(2) node(2)],'g','linewidth',2);
		else
			p(5) = plot([root(1) node(1)],[root(2) node(2)],'b','linewidth',2);
		end
		if k < size(combos,1)
			p(2) = plot(sumC,k,'m.','markersize',30);
		else
			p(3) = plot(sumC,k,'r.','markersize',30);
		end
		root = node;
	end
end
yticks(0:1:size(combos,1))
legend(p,{'Root','Branch','Leaf','$\xi^b_{i}=1$','$\xi^b_{i}=-1$'},'interpreter','latex','location','northeastoutside')
ylabel('i for $\xi^b_i$','interpreter','latex')
set(gca,'YDir','reverse','fontsize',14)
set(gca,'xtick',[])
set(gca,'xticklabel',[])

end