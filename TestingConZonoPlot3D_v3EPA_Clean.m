rng(3) % 2 was the working example
nG = 40; % 5 was the working example
nC = 15; % 1 was the working example
G = 2*rand(3,nG)-1;
c = zeros(3,1);
A = 2*rand(nC,nG)-1;
b = 2*rand(nC,1)-1;
obj = conZono(G,c,A,b);

% Still need to handle this case
% obj = conZono([eye(3) [1;1;1]],zeros(3,1),[1 1 1 0],[1]); 

%% MPT plot

figure;  hold on; axis square; view(3);
% tStart = tic;
% Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1),'He',[obj.A obj.b]);
% P = obj.c + obj.G*Box;
% plot(P,'color','g','alpha',0.1)
% % drawnow
% toc(tStart)

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
% plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok');view(3)
% Find opposite vertex (Need to add case to handle if save vertex is
% returned
dir = -[1 0 0];
[extreme] = getVertices(obj,model,params,dir);
[foundVerts,~] = addIfNew(foundVerts,extreme);
% plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok');view(3)
% Find 3rd vertex
nullDirs = null(foundVerts(1,:)-foundVerts(2,:))';
dir = nullDirs(1,:);
[extreme] = getVertices(obj,model,params,dir);
[foundVerts,~] = addIfNew(foundVerts,extreme);
% plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok');view(3)
% Find 4rd vertex
dir = nullDirs(2,:);
[extreme] = getVertices(obj,model,params,dir);
[foundVerts,~] = addIfNew(foundVerts,extreme);
% plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok');view(3)
%%%%%%%%%%%%%%% Need to handle case if this does not produce a 3-simplex
% Interior point
interiorPoint = sum(foundVerts(1:4,:))/4;
% plot3(interiorPoint(:,1),interiorPoint(:,2),interiorPoint(:,3),'or');view(3)
% Faces, vertices, and normals
[f] = nchoosek(1:4,3);
fNorms = zeros(0,3);
for i = 1:size(f,1)
    [fN,fC] = faceNormal(foundVerts(f(i,1:3),:),interiorPoint);
    fNorms(i,:)= fN;
end
% Expand polytope
maxIter = 1e6;
fIndex = 1;
nVerts = 4;
for i = 1:maxIter % Good up to 11
    [dir,~] = faceNormal(foundVerts(f(fIndex,1:3),:),interiorPoint);
    [extreme] = getVertices(obj,model,params,dir);
    [foundVerts,isNew] = addIfNew(foundVerts,extreme);
    if isNew
        nVerts = nVerts + 1;
%         plot3(extreme(:,1),extreme(:,2),extreme(:,3),'ok');view(3)
        % Find all faces that see the new vertex
        facesThatSee = findFacesThatSeeNewVertex(fNorms,foundVerts(f(:,1),:),extreme);
        sharedEdges = findSharedEdges(find(facesThatSee),f);
        facesToRemove = [];
        for indx = find(facesThatSee)
            removeCurrentFace = 0;
            if ismember(nVerts,f(indx,:))
                continue
            end
%             allCombos = nchoosek(f(indx,f(indx,:)~=0),2);
            nVertsRow = find(f(indx,:),1,'last');
            startVerts = f(indx,1:nVertsRow);
            endVerts = [f(indx,2:nVertsRow) f(indx,1)];
            allCombos = [startVerts' endVerts'];
%             nchoosek(f(indx,f(indx,:)~=0),2)
%             [startVerts' endVerts']
            allCombos = setdiff(sort(allCombos,2),sort(sharedEdges,2),'rows'); 
            newFaces = [allCombos nVerts*ones(size(allCombos,1),1)];
            minIndices = [];
            if isempty(newFaces)
                removeCurrentFace = 1;
            end
            for j = 1:size(newFaces,1)
                [fN,~] = faceNormal(foundVerts(newFaces(j,1:3),:),interiorPoint);
                [minVal,minIndx] = min(sum(abs(fNorms-fN),2));
                isNewFace = minVal>=1e-6;                   % Tolerance
                if isNewFace % Create new face
                    fNorms = [fNorms;fN];
                    f = [f;newFaces(j,:) zeros(1,size(f,2)-3)];
                    removeCurrentFace = 1;
                else % Combine with existing face
                    if ~ismember(nVerts,f(minIndx,:))
                        minIndices = [minIndices;minIndx];
                        n = find(f(minIndx,:),1,'last');
                        if n + 1 > size(f,2)
                            f = [f zeros(size(f,1),1)];
                        end
                        f(minIndx,n+1) = nVerts;
                        f(minIndx,:) = reOrderVerts(f(minIndx,:),foundVerts,interiorPoint);
                        removeCurrentFace = 1;
                    end
                end
            end
            if removeCurrentFace && ~ismember(indx,minIndices)
                facesToRemove = [facesToRemove; indx];
            end
        end
        f(facesToRemove,:) = [];
        fNorms(facesToRemove,:) = [];
    else
        fIndex = fIndex + 1;
    end
    if fIndex > size(f,1)
        break
    end
        
end
v = foundVerts;
f_ = f;
f_(f==0) = nan;
patch('Faces',f_,'Vertices',v,'FaceColor','m','FaceAlpha',1)
for m = 1:size(f_,1)
    patch('Faces',f_(m,:),'Vertices',v,'FaceColor','m','FaceAlpha',0.5)
end
% drawnow
plotbrowser
toc(tStart)

P
nfacesAndVerts = [size(f,1) nVerts]

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
allEdges = [];
for i = 1:length(faceIndices)
    nVerts = find(f(faceIndices(i),:),1,'last');
    startVerts = f(faceIndices(i),1:nVerts);
    endVerts = [f(faceIndices(i),2:nVerts) f(faceIndices(i),1)];
    allEdges = [allEdges; startVerts' endVerts'];
end
% startVerts = f(faceIndices,1:end);
% endVerts = [f(faceIndices,2:end) f(faceIndices,1)];
% allEdges = [reshape(startVerts,[],1) reshape(endVerts,[],1)];
allEdges = sort(allEdges,2);
[uniqueEdges,uniqueIndices] = unique(allEdges,'rows');
repeatedIndices = [1:size(allEdges,1)]';
repeatedIndices(uniqueIndices) = []; 
sharedEdges = allEdges(repeatedIndices,:);
end

function f = reOrderVerts(f,foundVerts,interiorPoint)
    nVerts = find(f(1,:),1,'last');
    [fN,fC] = faceNormal(foundVerts(f(1,1:3),:),interiorPoint);
    subSpaceVecs = foundVerts(f(1,1:nVerts),:)*null(fN);
    subCenter = fC*null(fN);
    angles = atan2(subSpaceVecs(:,1)-subCenter(:,1),subSpaceVecs(:,2)-subCenter(:,2));
    [~,order] = sort(angles);
    f(1,1:nVerts) = f(1,order);
end