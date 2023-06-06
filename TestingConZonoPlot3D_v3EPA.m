% rng(4) % Works 
% rng(1) % Does not work
rng(2)
nG = 4;
nC = 1;
G = 2*rand(3,nG)-1;
c = zeros(3,1);
A = 2*rand(nC,nG)-1;
b = 2*rand(nC,1)-1;
obj = conZono(G,c,A,b);

% obj = conZono([eye(3) [1;1;1]],zeros(3,1),[1 1 1 0],[1]);

%% CORA plot
% CORA = conZonotope(c,G,A,b);
% figure;
% tStart = tic;
% plot(CORA)
% drawnow
% toc(tStart)

% MPT plot

figure;  hold on; axis square
tStart = tic;
Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1),'He',[obj.A obj.b]);
P = obj.c + obj.G*Box;
plot(P,'color','g','alpha',0.1)
drawnow
toc(tStart)

%%
% figure; hold on; axis square
tStart = tic;
% Gurobi LP model setup
model.A = sparse(obj.A);
model.rhs = [obj.b];
model.ub = ones(obj.nG,1);
model.lb = -ones(obj.nG,1);
model.sense = '=';
model.modelsense = 'max';
params.Threads = 1;	% number of cores -> 0 => will use however many it wants
params.outputflag = 0;
% Initialization
foundVerts = zeros(0,3);
% Find a 3-simplex
% First vertex
dir = [1 0 0];
[extreme] = getVertices(obj,model,params,dir);
[foundVerts,~] = addIfNew(foundVerts,extreme);
plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok');view(3)
% Find opposite vertex (Need to add case to handle if save vertex is
% returned
dir = -[1 0 0];
[extreme] = getVertices(obj,model,params,dir);
[foundVerts,~] = addIfNew(foundVerts,extreme);
plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok');view(3)
% Find 3rd vertex
nullDirs = null(foundVerts(1,:)-foundVerts(2,:))';
dir = nullDirs(1,:);
[extreme] = getVertices(obj,model,params,dir);
[foundVerts,~] = addIfNew(foundVerts,extreme);
plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok');view(3)
% Find 4rd vertex
dir = nullDirs(2,:);
[extreme] = getVertices(obj,model,params,dir);
[foundVerts,~] = addIfNew(foundVerts,extreme);
plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok');view(3)
%%%%%%%%%%%%%%% Need to handle case if this does not produce a 3-simplex
% Interior point
interiorPoint = sum(foundVerts(1:4,:))/4;
plot3(interiorPoint(:,1),interiorPoint(:,2),interiorPoint(:,3),'or');view(3)
% Faces, vertices, and normals
[f] = nchoosek(1:4,3);
v = foundVerts;
fNorms = zeros(0,3);
for i = 1:size(f,1)
    patch('Faces',f(i,:),'Vertices',v,'FaceColor','m','FaceAlpha',0.1)
    [fN,fC] = faceNormal(foundVerts(f(i,1:3),:),interiorPoint);
    fNorms(i,:)= fN;
    quiver3(fC(:,1),fC(:,2),fC(:,3), ...
                fN(:,1),fN(:,2),fN(:,3),1.5,'color','b');
end
% Expand polytope
maxIter = 1e6;
fIndex = 1;
nVerts = 4;
for i = 1:2%maxIter
    [dir,~] = faceNormal(foundVerts(f(fIndex,1:3),:),interiorPoint);
    [extreme] = getVertices(obj,model,params,dir);
    [foundVerts,isNew] = addIfNew(foundVerts,extreme);
    if isNew
        nVerts = nVerts + 1;
        plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok');view(3)
        % Find all faces that see the new vertex
        facesThatSee = findFacesThatSeeNewVertex(fNorms,foundVerts(f(:,1),:),extreme);
        sharedEdges = findSharedEdges(find(facesThatSee),f);
        removeCurrentFace = 0;
        for indx = find(facesThatSee)
%             removeCurrentFace = 0;
            allCombos = nchoosek(f(fIndex,f(fIndex,:)~=0),2);
            allCombos = setdiff(allCombos,sharedEdges,'rows'); 
            newFaces = [allCombos nVerts*ones(size(allCombos,1),1)];
            f
            for j = 1:size(newFaces,1)
                [fN,~] = faceNormal(foundVerts(newFaces(j,1:3),:),interiorPoint);
                [minVal,minIndx] = min(sum(abs(fNorms-fN),2));
                isNewFace = minVal>=1e-6;                   % Tolerance
                if isNewFace % Create new face
                    fNorms = [fNorms;fN];
                    f = [f;newFaces(j,:) zeros(1,size(f,2)-3)];
%                     removeCurrentFace = 1;
                else % Combine with existing face
                    if ~ismember(nVerts,f(minIndx,:))
                        n = find(f(minIndx,:),1,'last');
                        if n + 1 > size(f,2)
                            f = [f zeros(size(f,1),1)];
                        end
                        f(minIndx,n+1) = nVerts;
                        f(minIndx,:) = reOrderVerts(f(minIndx,:),foundVerts,interiorPoint);
                    end
                end
                removeCurrentFace = 1;
            end
%             if removeCurrentFace
%                 f(indx,:) = [];
%                 fNorms(indx,:) = [];
%             end
        end
%         f
% find(facesThatSee)
%         for indx = find(facesThatSee)
%             
% 
% 
%             if indx == fIndex % Create new faces
%                 allCombos = nchoosek(f(fIndex,f(fIndex,:)~=0),2);
%                 newFaces = [allCombos nVerts*ones(size(allCombos,1),1)];
%                 for j = 1:size(newFaces,1)
%                     f = [f;newFaces(j,:) zeros(1,size(f,2)-3)];
%                     [fN,~] = faceNormal(foundVerts(f(j,1:3),:),interiorPoint);
%                     fNorms = [fNorms;fN]; % Need to check if new normal
%                 end
%                 removeCurrentFace = 1;
%             else % Add to existing faces
%                 n = find(f(indx,:),1,'last');
%                 if n + 1 > size(f,2)
%                     f = [f zeros(size(f,1),1)];
%                 end
%                 f(indx,n+1) = nVerts;
%             end
%         end
        if removeCurrentFace
            f(fIndex,:) = [];
            fNorms(fIndex,:) = [];
        end
    else
        fIndex = fIndex + 1;
    end
        
end

% Reorder vertices in faces
% for i = 1:size(f,1)
%     f(i,:) = reOrderVerts(f(i,:),foundVerts,interiorPoint);
% end
v = foundVerts;
f_ = f;
f_(f==0) = nan;
% patch('Faces',f_,'Vertices',v,'FaceColor','m','FaceAlpha',1)
for m = 1:size(f_,1)
    patch('Faces',f_(m,:),'Vertices',v,'FaceColor','m','FaceAlpha',0.5)
end
drawnow
toc(tStart)

%%
% % figure; hold on; axis square
% tStart = tic;
% nLeaves = 1;
% % Gurobi LP model setup
% model.A = sparse(obj.A);
% model.rhs = [obj.b];
% model.ub = ones(obj.nG,1);
% model.lb = -ones(obj.nG,1);
% model.sense = '=';
% model.modelsense = 'max';
% % define the solver options to look for all feasible solutions
% params.Threads = 1;	% number of cores -> 0 => will use however many it wants
% params.outputflag = 0;
% % Initialize with first 2 vertices
% maxVerts = 100;
% foundVerts = zeros(maxVerts*nLeaves,3);
% % vertOrder = zeros(maxVerts*nLeaves,1);
% % treeStruct = zeros(maxVerts*nLeaves,3);
% searchedDirs = zeros(0,3);
% % totalVerts = 0;
% foundVertDirs = zeros(maxVerts*nLeaves,2);
% firstIndices = [1:maxVerts:maxVerts*nLeaves]';
% % % Find first vertex
% dir = [1 0 0];
% searchedDirs = [searchedDirs;normalize(dir,'norm')];
% [extreme] = getVertices(obj,model,params,dir);
% foundVerts(firstIndices,:) = extreme;
% [a,e,~] = cart2sph(extreme(:,1),extreme(:,2),extreme(:,3));
% foundVertDirs(firstIndices,:) = [a e];
% % plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok')
% view(3)
% % Find opposite vertex
% dir = -[1 0 0];
% searchedDirs = [searchedDirs;normalize(dir,'norm')];
% [extreme] = getVertices(obj,model,params,dir);
% foundVerts(firstIndices+1,:) = extreme;
% [a,e,~] = cart2sph(extreme(:,1),extreme(:,2),extreme(:,3));
% foundVertDirs(firstIndices+1,:) = [a e];
% % plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok')
% view(3)
% % Find 3rd vertex
% nullDirs = null(foundVerts(1,:)-foundVerts(2,:))';
% % dir = [0 -1 0];
% dir = nullDirs(1,:);
% searchedDirs = [searchedDirs;normalize(dir,'norm')];
% [extreme] = getVertices(obj,model,params,dir);
% foundVerts(firstIndices+2,:) = extreme;
% [a,e,~] = cart2sph(extreme(:,1),extreme(:,2),extreme(:,3));
% foundVertDirs(firstIndices+2,:) = [a e];
% % plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok')
% view(3)
% % Interior point
% interiorPoint = sum(foundVerts(1:3,:))/3;
% % plot3(interiorPoint(:,1),interiorPoint(:,2),interiorPoint(:,3),'og')
% view(3)
% nVerts = 3;
% % Initial two faces
% f = [1 2 3; 2 1 3];
% fNorms = [];
% [fNorms(1,:),~] = faceNormal(foundVerts(f(1,1:3),:),interiorPoint);
% [fNorms(2,:),~] = faceNormal(foundVerts(f(2,1:3),:),interiorPoint);
% % f = [1 2 3];
% maxIter = 1e6;
% fIndex = 1;
% for i = 1:100%maxIter
%     [fN,~] = faceNormal(foundVerts(f(fIndex,1:3),:),interiorPoint);
%     dir = fN;
%     searchedDirs = [searchedDirs;normalize(dir,'norm')];
%     [extreme] = getVertices(obj,model,params,dir);
%     [a,e,~] = cart2sph(extreme(:,1),extreme(:,2),extreme(:,3));
%     isNew = min(sum(abs(foundVertDirs-[a e]),2))>=1e-6;                   % Tolerance
%     if ~isNew
%         [~,vertIndx] = min(sum(abs(foundVertDirs-[a e]),2));
%         vertInFace = min(abs(f(fIndex,:)-vertIndx))<=1e-6;
%         if ~vertInFace
%             f
%             fIndex
%             allCombos = nchoosek(f(fIndex,f(fIndex,:)~=0),2);
%             newFaces = [allCombos vertIndx*ones(size(allCombos,1),1)];
%             moddedFaces = 0;
%             for j = 1:size(newFaces,1)
%                 [fNorm,fCenter] = faceNormal(foundVerts(newFaces(j,:),:),interiorPoint);
%                 [minVal,minIndx] = min(sum(abs(fNorms-fNorm),2));
%                 isNewFace = minVal>=1e-6;                   % Tolerance
%                 if isNewFace % Create new face
%                     moddedFaces = 1;
% %                     f(fIndex,:) = [];
% %                     fNorms(fIndex,:) = [];
%                     fNorms = [fNorms;fNorm];
%                     f = [f;newFaces(j,:) zeros(1,size(f,2)-3)];
%                 else % Combine with existing face
% % %                     %                 vertInFace = min(abs(f(minIndx,:)-vertIndx))<=1e-6;
% % %                     % %                 if ~vertInFace % Add to face
%                     n = find(f(minIndx,:),1,'last');
%                     if n == size(f,2)
%                         f = [f zeros(size(f,1),1)];
%                     end
%                     if ~isempty(setdiff(newFaces(j,:),f(minIndx,:)))
%                         moddedFaces = 1;
%                         f(minIndx,n+1) = setdiff(newFaces(j,:),f(minIndx,:));
%                         % Order vertices
%                         subSpaceVecs = foundVerts(f(minIndx,1:n+1),:)*null(fNorm);
%                         subCenter = fCenter*null(fNorm);
%                         angles = atan2(subSpaceVecs(:,1)-subCenter(:,1),subSpaceVecs(:,2)-subCenter(:,2));
%                         [~,order] = sort(angles);
%                         f(minIndx,1:n+1) = f(minIndx,order);
%                     end
%                 end
%                 if moddedFaces
%                     f(fIndex,:) = [];
%                     fNorms(fIndex,:) = [];
%                 end
%             end
%         else
%             fIndex = fIndex + 1;
%         end
%             
% %         end
%         if fIndex == size(f,1)
%             break
% %         elseif ~addedNewFaces
% %             fIndex = fIndex + 1;
%         end
%     else
%         nVerts = nVerts + 1;
%         foundVerts(firstIndices+nVerts-1,:) = extreme;
%         foundVertDirs(firstIndices+nVerts-1,:) = [a e];
% %         plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok')
%         view(3)
%         allCombos = nchoosek(f(fIndex,f(fIndex,:)~=0),2);
%         newFaces = [allCombos nVerts*ones(size(allCombos,1),1)];
%         f(fIndex,:) = [];
%         fNorms(fIndex,:) = [];
%         for j = 1:size(newFaces,1)
%             [fNorm,fCenter] = faceNormal(foundVerts(newFaces(j,:),:),interiorPoint);
%             [minVal,minIndx] = min(sum(abs(fNorms-fNorm),2));
%             isNewFace = minVal>=1e-6;                   % Tolerance
%             if isNewFace % Create new face
%                 fNorms = [fNorms;fNorm];
%                 f = [f;newFaces(j,:) zeros(1,size(f,2)-3)];
%             else % Combine with existing face
%                 n = find(f(minIndx,:),1,'last');
%                 verts2Add = setdiff(newFaces(j,:),f(minIndx,:));
%                 if n + length(verts2Add) > size(f,2)
%                     f = [f zeros(size(f,1),length(verts2Add))];
%                 end
%                 if ~isempty(verts2Add)
%                     f(minIndx,n+(1:length(verts2Add))) = verts2Add;
%                     % Order vertices
%                     subSpaceVecs = foundVerts(f(minIndx,1:n+1),:)*null(fNorm);
%                     subCenter = fCenter*null(fNorm);
%                     angles = atan2(subSpaceVecs(:,1)-subCenter(:,1),subSpaceVecs(:,2)-subCenter(:,2));
%                     [~,order] = sort(angles);
%                     f(minIndx,1:n+1) = f(minIndx,order);
%                 end
%             end
%             
%         end
% 
% %         for k = 1:size(f,1)
% %             [fN,fC] = faceNormal(foundVerts(f(k,1:3),:),interiorPoint);
% %             quiver3(fC(:,1),fC(:,2),fC(:,3), ...
% %                 fN(:,1),fN(:,2),fN(:,3),1.5,'color','b');
% %         end
% %         v = foundVerts;
% %         f_ = f;
% %         f_(f==0) = nan;
% %         for m = 1:size(f_,1)
% %             patch('Faces',f_(m,:),'Vertices',v,'FaceColor','m','FaceAlpha',0.1)
% %         end
% % %         patch('Faces',f_,'Vertices',v,'FaceColor','m','FaceAlpha',0.1)
% %         drawnow
% %         f
% %         fNorms
%     end
% end
% 
% 
% % for i = 1:size(f,1)
% %     [fN,fC] = faceNormal(foundVerts(f(i,:),:),interiorPoint);
% %     quiver3(fC(:,1),fC(:,2),fC(:,3), ...
% %         fN(:,1),fN(:,2),fN(:,3),1.5,'color','b');
% % end
% v = foundVerts;
% f_ = f;
% f_(f==0) = nan;
% patch('Faces',f_,'Vertices',v,'FaceColor','m','FaceAlpha',1)
% % for m = 1:size(f_,1)
% %     patch('Faces',f_(m,:),'Vertices',v,'FaceColor','m','FaceAlpha',1)
% % end
% drawnow
% toc(tStart)
% P
% nfacesAndVerts = [size(f,1) nVerts]
% 
% %%
% % m = 5;
% % [fN,~] = faceNormal(foundVerts(f(5,1:3),:),interiorPoint);
% % dir = fN;
% % [extreme] = getVertices(obj,model,params,dir);
% % [a,e,~] = cart2sph(extreme(:,1),extreme(:,2),extreme(:,3));
% % isNew = min(sum(abs(foundVertDirs-[a e]),2))>=1e-6;
% % [~,vertIndx] = min(sum(abs(foundVertDirs-[a e]),2));
% % vertInFace = min(abs(f(fIndex,:)-vertIndx))<=1e-6;
% % %%
% % m = 5
% % [fNorm,fCenter] = faceNormal(foundVerts(f(m,1:3),:),interiorPoint)f
% % subSpaceVecs = foundVerts(f(m,1:3),:)*null(fNorm)
% % subCenter = fCenter*null(fNorm)
% % angles = atan2(subSpaceVecs(:,1)-subCenter(:,1),subSpaceVecs(:,2)-subCenter(:,2))
% % [~,order] = sort(angles)
% % figure; hold on; axis square
% % patch('Faces',f_(m,order),'Vertices',v,'FaceColor','m','FaceAlpha',1)
% % [fN,fC] = faceNormal(foundVerts(f(m,1:3),:),interiorPoint);
% % plot3(fC(:,1),fC(:,2),fC(:,3),'og')
% % quiver3(fC(:,1),fC(:,2),fC(:,3), ...
% %     fN(:,1),fN(:,2),fN(:,3),1.5,'color','b');
%%
function [extreme] = getVertices(obj,model,params,dir)
model.obj = dir*obj.G;
result = gurobi(model,params);
extreme = [obj.G*result.x + obj.c]';
end

function [foundVerts,isNew] = addIfNew(foundVerts,extreme)
    isNew = 0;
    if isempty(foundVerts)
        foundVerts = [foundVerts;extreme];
        isNew = 1;
        return
    end
    minDiff = min(sum(abs(foundVerts-extreme),2));
    if minDiff >= 1e-6 % Threshold
        isNew = 1;
        foundVerts = [foundVerts;extreme];
    end
end

function [fN,fC] = faceNormal(verts,interiorPoint)
vecA = verts(2,:) - verts(1,:);
vecB = verts(3,:) - verts(1,:);

fN = cross(vecA,vecB);
fC = sum(verts)/size(verts,1);
dotProd = fN*[fC-interiorPoint]';
if dotProd < 0 % Face normals should always point outward
    fN = -fN;
end
fN = normalize(fN,'norm');
end

function facesThatSee = findFacesThatSeeNewVertex(fNorms,vertOnFace,extreme)
    distances = dot((extreme - vertOnFace)',fNorms');
    facesThatSee = distances >= -1e-6; % Threshold
end

function sharedEdges = findSharedEdges(faceIndices,f)
startVerts = f(faceIndices,1:end-1);
endVerts = f(faceIndices,2:end);
allEdges = [reshape(startVerts,[],1) reshape(endVerts,[],1)];
[uniqueEdges,uniqueIndices] = unique(allEdges,'rows');
repeatedIndices = [1:size(allEdges,1)]';
repeatedIndices(uniqueIndices) = []; 
sharedEdges = allEdges(repeatedIndices,:);
end

function f = reOrderVerts(f,foundVerts,interiorPoint)
    for j = 1:find(f(1,:),1,'last')-2
        vertA = foundVerts(f(1,j),:);
        vertB = foundVerts(f(1,j+1),:);
        vertC = foundVerts(f(1,j+2),:);
        vecA = vertB - vertA;
        vecB = vertC - vertB;
        fN = cross(vecA,vecB);
        fC = (vertA + vertB + vertC)/3;
        dotProd = fN*[fC-interiorPoint]';
        if dotProd > 0 % Switch order of vertex in face
            f1 = f(1,j+1);
            f2 = f(1,j+2);
            f(1,j+1) = f2;
            f(1,j+2) = f1;
        end
    end
end