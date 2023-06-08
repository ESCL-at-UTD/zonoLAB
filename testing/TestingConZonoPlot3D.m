rng(4)
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
nLeaves = 1;
% Gurobi LP model setup
model.A = sparse(obj.A);
model.rhs = [obj.b];
model.ub = ones(obj.nG,1);
model.lb = -ones(obj.nG,1);
model.sense = '=';
model.modelsense = 'max';
% define the solver options to look for all feasible solutions
params.Threads = 1;	% number of cores -> 0 => will use however many it wants
params.outputflag = 0;
% Initialize with first 2 vertices
maxVerts = 100;
foundVerts = zeros(maxVerts*nLeaves,3);
% vertOrder = zeros(maxVerts*nLeaves,1);
% treeStruct = zeros(maxVerts*nLeaves,3);
searchedDirs = zeros(0,3);
% totalVerts = 0;
foundVertDirs = zeros(maxVerts*nLeaves,2);
firstIndices = [1:maxVerts:maxVerts*nLeaves]';
% % Find first vertex
dir = [1 0 0];
searchedDirs = [searchedDirs;normalize(dir,'norm')];
[extreme] = getVertices(obj,model,params,dir);
foundVerts(firstIndices,:) = extreme;
[a,e,~] = cart2sph(extreme(:,1),extreme(:,2),extreme(:,3));
foundVertDirs(firstIndices,:) = [a e];
plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok')
view(3)
% Find opposite vertex
dir = -[1 0 0];
searchedDirs = [searchedDirs;normalize(dir,'norm')];
[extreme] = getVertices(obj,model,params,dir);
foundVerts(firstIndices+1,:) = extreme;
[a,e,~] = cart2sph(extreme(:,1),extreme(:,2),extreme(:,3));
foundVertDirs(firstIndices+1,:) = [a e];
plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok')
view(3)
% Find 3rd vertex
nullDirs = null(foundVerts(1,:)-foundVerts(2,:))';
% dir = [0 -1 0];
dir = nullDirs(1,:);
searchedDirs = [searchedDirs;normalize(dir,'norm')];
[extreme] = getVertices(obj,model,params,dir);
foundVerts(firstIndices+2,:) = extreme;
[a,e,~] = cart2sph(extreme(:,1),extreme(:,2),extreme(:,3));
foundVertDirs(firstIndices+2,:) = [a e];
plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok')
view(3)
% Interior point
interiorPoint = sum(foundVerts(1:3,:))/3;
plot3(interiorPoint(:,1),interiorPoint(:,2),interiorPoint(:,3),'og')
view(3)
nVerts = 3;
% Initial two faces
f = [1 2 3; 2 1 3];
fNorms = [];
[fNorms(1,:),~] = faceNormal(foundVerts(f(1,1:3),:),interiorPoint);
[fNorms(2,:),~] = faceNormal(foundVerts(f(2,1:3),:),interiorPoint);
% f = [1 2 3];
maxIter = 1e6;
fIndex = 1;
for i = 1:100%maxIter
    [fN,~] = faceNormal(foundVerts(f(fIndex,1:3),:),interiorPoint);
    dir = fN;
    searchedDirs = [searchedDirs;normalize(dir,'norm')];
    [extreme] = getVertices(obj,model,params,dir);
    [a,e,~] = cart2sph(extreme(:,1),extreme(:,2),extreme(:,3));
    isNew = min(sum(abs(foundVertDirs-[a e]),2))>=1e-6;                   % Tolerance
    if ~isNew
        if fIndex == size(f,1)
            break
        else
            fIndex = fIndex + 1;
        end
    else
        nVerts = nVerts + 1;
        foundVerts(firstIndices+nVerts-1,:) = extreme;
        foundVertDirs(firstIndices+nVerts-1,:) = [a e];
        plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok')
        view(3)
        allCombos = nchoosek(f(fIndex,f(fIndex,:)~=0),2);
        newFaces = [allCombos nVerts*ones(size(allCombos,1),1)];
        f(fIndex,:) = [];
        fNorms(fIndex,:) = [];
        for j = 1:size(newFaces,1)
            [fNorm,~] = faceNormal(foundVerts(newFaces(j,:),:),interiorPoint);
            [minVal,minIndx] = min(sum(abs(fNorms-fNorm),2));
            isNewFace = minVal>=1e-6;                   % Tolerance
            if isNewFace % Create new face
                fNorms = [fNorms;fNorm];
                f = [f;newFaces(j,:) zeros(1,size(f,2)-3)];
            else % Combine with existing face
                n = find(f(minIndx,:),1,'last');
                if n == size(f,2)
                    f = [f zeros(size(f,1),1)];
                end
                if ~isempty(setdiff(newFaces(j,:),f(minIndx,:)))
                    f(minIndx,n+1) = setdiff(newFaces(j,:),f(minIndx,:));
                    % Order vertices
                    subSpaceVecs = foundVerts(f(minIndx,1:n+1),:)*null(fNorm);
                    angles = atan2(subSpaceVecs(:,1),subSpaceVecs(:,2));
                    [~,order] = sort(angles);
                    f(minIndx,1:n+1) = f(minIndx,order);
                end
            end
            
        end


%         for k = 1:size(f,1)
%             [fN,fC] = faceNormal(foundVerts(f(k,1:3),:),interiorPoint);
%             quiver3(fC(:,1),fC(:,2),fC(:,3), ...
%                 fN(:,1),fN(:,2),fN(:,3),1.5,'color','b');
%         end
%         v = foundVerts;
%         f_ = f;
%         f_(f==0) = nan;
%         patch('Faces',f_,'Vertices',v,'FaceColor','m','FaceAlpha',0.1)
%         drawnow
    end
end


% for i = 1:size(f,1)
%     [fN,fC] = faceNormal(foundVerts(f(i,:),:),interiorPoint);
%     quiver3(fC(:,1),fC(:,2),fC(:,3), ...
%         fN(:,1),fN(:,2),fN(:,3),1.5,'color','b');
% end
v = foundVerts;
f_ = f;
f_(f==0) = nan;
% patch('Faces',f_,'Vertices',v,'FaceColor','m','FaceAlpha',1)
for m = 1:size(f_,1)
    patch('Faces',f_(m,:),'Vertices',v,'FaceColor','m','FaceAlpha',1)
end
drawnow
toc(tStart)
% % Compute vertex angles
% vertDiffs = firstVert - centerVert;
% referenceAngles = atan2(vertDiffs(:,2),vertDiffs(:,1));
% foundVertDirs(firstIndices) = zeros(nLeaves,2);
% % Farthest CW vert
% foundVerts(firstIndices+1,:) = endVert;
% foundVertDirs(firstIndices+1) = -pi;
% treeStruct(firstIndices,3) = firstIndices+1;
% treeStruct(firstIndices+1,:) = [firstIndices zeros(nLeaves,2)];
% % Farthest CW vert
% foundVerts(firstIndices+2,:) = endVert;
% foundVertDirs(firstIndices+2) = pi;
% treeStruct(firstIndices,2) = firstIndices+2;
% treeStruct(firstIndices+2,:) = [firstIndices zeros(nLeaves,2)];
% % Find the remaining vertices
% nVerts = 3*ones(nLeaves,1);
% for i = 1:nLeaves
%     [vertOrder(firstIndices(i)+[0:nVerts(i)-1]),~] = listInOrder(treeStruct,firstIndices(i),zeros(nVerts(i),1),1);
% end
% maxIter = 1e6;
% indx = 1;
% j = 1;
% for i = 1:maxIter
%     vertA = foundVerts(vertOrder(indx),:);
%     vertB = foundVerts(vertOrder(indx+1),:);
%     vertDiff = vertA - vertB;
%     dir = -[-vertDiff(2) vertDiff(1)];
%     isNewdir = min(sum(abs(normalize(dir,'norm')-searchedDirs),2))>=1e-6;   % Tolerance
%     if isNewdir
%         searchedDirs = [searchedDirs;normalize(dir,'norm')];                % Growing list
%         [extreme] = getVertices(obj,model,params,dir);
%         vertDiffs = extreme - centerVert;
%         vertAngles = atan2(vertDiffs(:,2),vertDiffs(:,1));
%         newVertDirs = vertAngles-referenceAngles;
%         newVertDirs(newVertDirs<-pi) = newVertDirs(newVertDirs<-pi) + 2*pi;
%         newVertDirs(newVertDirs>pi) = newVertDirs(newVertDirs>pi) - 2*pi;
%         orderDirs = foundVertDirs(vertOrder(firstIndices(j)+[0:nVerts(j)-1]));
%         isNew = min(abs(newVertDirs(j)-orderDirs))>=1e-6;                   % Tolerance
%     else
%         isNew = 0;
%     end
% 
%     if isNew
%         nVerts(j) = nVerts(j) + 1;
%         foundVerts(firstIndices(j)+nVerts(j)-1,:) = extreme(j,:);
%         foundVertDirs(firstIndices(j)+nVerts(j)-1) = newVertDirs(j);
%         if treeStruct(vertOrder(indx),3) == 0
%             treeStruct(vertOrder(indx),3) = firstIndices(j)+nVerts(j)-1;
%             treeStruct(firstIndices(j)+nVerts(j)-1,:) = [vertOrder(indx) 0 0];
%         else
%             treeStruct(vertOrder(indx+1),2) = firstIndices(j)+nVerts(j)-1;
%             treeStruct(firstIndices(j)+nVerts(j)-1,:) = [vertOrder(indx+1) 0 0];
%         end
%         [vertOrder(firstIndices(j)+[0:nVerts(j)-1]),~] = listInOrder(treeStruct,firstIndices(j),zeros(nVerts(j),1),1);
% %         plot(extreme(j,1),extreme(j,2),'ok')
%     else
%         if vertOrder(indx+1) == firstIndices(j)+1
%             totalVerts(j) = sum(nVerts);
%             if j == nLeaves
%                 break
%             else
%                 j = j + 1;
%                 indx = firstIndices(j);
%             end
%         else
%             indx = indx + 1;
%         end
%     end
% end
% nVerts = nVerts - 1; % Added to remove extra vertex (May not work with hybZono)
% v = foundVerts(vertOrder(vertOrder~=0),:);
% f = nan*ones(nLeaves,max(nVerts));
% v = v(1:end-1,:); % Added to remove extra vertex (May not work with hybZono)
% % f = f(:,1:end-1); % Added to remove extra vertex (May not work with hybZono)
% count = 1;
% for i2 = 1:nLeaves
%     f(i2,1:nVerts(i2)) = count+[0:nVerts(i2)-1];
%     count = count + nVerts(i2);
% end
% v = foundVerts;
% patch('Faces',f,'Vertices',v,'FaceColor','m','FaceAlpha',1)
% drawnow
% toc(tStart)
% 
%%
function [extreme] = getVertices(obj,model,params,dir)
model.obj = dir*obj.G;
result = gurobi(model,params);
extreme = [obj.G*result.x + obj.c]';
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
% 
% function [vertOrder,indx] = listInOrder(treeStruct,row,vertOrder,indx)
%     if treeStruct(row,2) ~= 0
%         [vertOrder,indx] = listInOrder(treeStruct,treeStruct(row,2),vertOrder,indx);
%     end
%     vertOrder(indx) = row;
%     indx = indx + 1;
%     if treeStruct(row,3) ~= 0
%         [vertOrder,indx] = listInOrder(treeStruct,treeStruct(row,3),vertOrder,indx);
%     end   
% end