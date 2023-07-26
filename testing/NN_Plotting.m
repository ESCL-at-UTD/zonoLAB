% load HybZono_sin_cos_2_20_10_10_1_matrices.mat
% Z = hybZono(Gc,Gb,c,Ac,Ab,b);
% load HybZono_sin_cos_2_20_10_10_1_Plot_Data.mat
% load HybZono_sin_cos_2_20_20_20_1.mat
load HybZono_sin_cos_2_10_10_1.mat
Z = NN_HybZono;


optSolver = solverOptions('milpSolver','gurobi');
optSolver.lpSolver = 'gurobi';
figure; hold on
tStart = tic;
opt = plotOptions('FaceColor','r','FaceAlpha',1,'EdgeAlpha',1,'SolverOpts',optSolver);
opt.Display = 'individual';
[v,f] = plot(Z,opt);
drawnow
toc(tStart)
%%
faceNormals = zeros(size(f,1),3);
facePoints = zeros(size(f,1),3);
for i = 1:size(f,1)
    verts = v(f(i,1:3),:);
    vecA = verts(2,:) - verts(1,:);
    vecB = verts(3,:) - verts(1,:);

    fN = cross(vecA,vecB);
    if fN(3) < 0 % Want face normals to point upward
        fN = -fN;
    end
    faceNormals(i,:) = normalize(fN,'norm');
    facePoints(i,:) = sum(verts)/size(verts,1);
end
quiver3(facePoints(:,1),facePoints(:,2),facePoints(:,3),...
            faceNormals(:,1),faceNormals(:,2),faceNormals(:,3),'k')

%%
fKeep(:,1) = 1:size(f,1);
for i = 1:size(f,1)-1
    distances = [];
    for j = i+1:size(f,1)
        distances(j) = dot((facePoints(j,:) - facePoints(i,:))',faceNormals(i,:)');
    end
    isNegative = max(distances < 0);
    if isNegative
        fKeep(i) = 0;
    end
end
% relations = faceNormals(:,:)*facePoints(:,:)';
% 
% greater = relations > repmat(diag(relations),1,size(f,1));
% fKeep(:,1) = 1:size(f,1);
% 
% for i = 18:size(f,1)
%     if fKeep(i) == 0
%         continue
%     end
%     i
% %     indx = faceNormals(i,:)*facePoints(:,:)' > faceNormals(i,:)*facePoints(i,:)';
% %     fKeep(find(indx)) = 0;
%     for j = i+1:size(f,1)
%         if (greater(i,j) == 0) %&& (greater(j,i) == 1)
%             fKeep(j) = 0;
%         end
%     end
% 
% end

%%
figure;hold on
for indx = 1:size(f,1)
    if fKeep(indx) == 0
        continue
    end
    P.Faces = f(indx,:);
    P.Vertices = v;
    patch(P,'FaceAlpha',0.5)
end
view(3)

%%
for i = setdiff(1:size(v,1),f(139,:))
distances(i) = dot((v(i,:) - facePoints(139,:))',faceNormals(139,:)')
end
min(distances)