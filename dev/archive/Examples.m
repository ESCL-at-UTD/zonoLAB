%% Test the construction of hybrid zonotopes with different numbers of inputs.
c = zeros(2,1);
Gc = [eye(2) ones(2,1)];
Gb = eye(2);
Ac = ones(1,3);
Ab = ones(1,2);
b = 1;
Z1 = hybZono(c)
Z2 = hybZono(Gc,c)
Z3 = hybZono(Gc,Gb,c)
Z4 = hybZono(Gc,c,Ac,b)
Z6 = hybZono(Gc,Gb,c,Ac,Ab,b)
% Z5 = hybZono(Gc,Gb,c,Ac,b) % Error demonstration
%% Test the display of multiple hybrid zonotopes in a table.
Z_all(1) = Z1;
Z_all(2) = Z2;
Z_all(3) = Z3;
Z_all(4) = Z4;
Z_all(5) = Z6;
Z_all
%% Identify leaves and binary tree and plot binary tree
Z6.getLeaves
Z6.getBinTree;
figure;
Z6.plotBinaryTree
%% Plot hybrid zonotopes
figure;
plot_MPT(Z6,'r',0.5)
%% Check if set is empty
Z6.isempty
Z7 = hybZono(Gc,Gb,c,Ac,Ab,-10);
Z7.isempty
%% Check for non-empty intersection of hybrid zonotope and halfspace
h = [1 2];
f = [1];
H = Polyhedron('A',h,'b',f,'lb',[-5 -5],'ub',[5 5]);
figure;hold on
plot(H)
plot_MPT(Z6,'r',0.1)
Z8 = Z6;
Z8.c = -[1;1];
plot_MPT(Z8,'b',0.1)
checkZintH(Z6,h,f)
checkZintH(Z8,h,f)



