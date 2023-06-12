% load HybZono_sin_cos_2_20_10_10_1_matrices.mat
% Z = hybZono(Gc,Gb,c,Ac,Ab,b);
% load HybZono_sin_cos_2_20_10_10_1_Plot_Data.mat
load HybZono_sin_cos_2_20_20_20_1.mat
Z = NN_HybZono;


optSolver = solverOptions('milpSolver','gurobi');
optSolver.lpSolver = 'gurobi';
figure;
tStart = tic;
opt = plotOptions('FaceColor','r','FaceAlpha',1,'EdgeAlpha',1,'SolverOpts',optSolver);
[v,f] = plot(Z,opt);
drawnow
toc(tStart)

% %%
% H = Z;
% combos = getLeaves(H,optSolver);
% nLeaves = size(combos,2);
% v = [];
% f = nan(nLeaves,1);
% figure;hold on
% for i = 1:nLeaves
%     i
%     clear P
%     tStart = tic;
%     obj = conZono(H.Gc,H.c+H.Gb*combos(:,i),H.Ac,H.b-H.Ab*combos(:,i));
%     T = null(obj.A,'r');
%     s = pinv(obj.A)*obj.b;
%     Box = Polyhedron('H',[T ones(size(s,1),1)-s; -T ones(size(s,1),1)+s]);
%     
%     obj.c = obj.c+obj.G*s;
%     obj.G = obj.G*T;
%     indx = find(sum(abs(obj.G))>=1e-6);
%     obj.G = obj.G(:,indx);
%     Box = projection(Box,indx);
%     
%     P = obj.c + obj.G*Box;
%         
%     nV = size(P.V(:,1),1);
%     if nV > size(f,2)
%         f = [f nan(size(f,1),nV-size(f,2))];
%     end
%     dir = atan2((P.V(:,1)-mean(P.V(:,1))),(P.V(:,2)-mean(P.V(:,2))));
%     [~,indx] = sort(dir,'ascend');
%     nV_all = size(v,1);
%     v = [v;P.V(indx,:)];
%     f(i,1:nV) = nV_all+[1:nV];
%     toc(tStart)
% end
% tic
% patch('Faces',f,'Vertices',v,'FaceColor','red');
% xlabel('x');ylabel('y');zlabel('z')
% view(3)
% toc