% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return vertices and faces for a constrained zonotope in 3D
%   Syntax:
%       [v,f] = plotConZono3D(Z,optSolver)
%   Inputs:
%       Z - 3D constrained zonotope in CG-Rep (conZono object)
%       optSolver - solver options needed for linear propgram
%   Outputs:
%       v - nV x 3 matrix, each row denoting the x (first column), y (second column), 
%                          and z (third column) positions of the nV vertices
%       f - nF x nMax matrix, each row denoting the vertices (up to nMax) contained
%                          in the nF faces (padded with NaN if face
%                          contains less than nMax vertices)
%   Notes:
%       Not intended to be called directly by user.
%       Use [v,f] = plot(obj,varargin) instead (method of abstractZono)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [v,f] = plotConZono3D(obj,optSolver)

% Problem data for linear program (LP)
Aeq = sparse(obj.A);
beq = [obj.b];
lb = -ones(obj.nG,1);
ub =  ones(obj.nG,1);

try
    % Find first vertex
    dir = [1 0 0];
    searchedDirs(1,:) = normalize(dir,'norm');
    %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
    x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
    foundVerts(1,:) = [obj.G*x + obj.c]';

    % Find second (opposite) vertex
    dir = [-1 0 0];
    searchedDirs(2,:) = normalize(dir,'norm');
    %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
    x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
    foundVerts(2,:) = [obj.G*x + obj.c]';

    if (foundVerts(1,:)-foundVerts(2,:) <= 1e-12) % Same vertex
        % Find first vertex
        dir = [0 1 0];
        searchedDirs(1,:) = normalize(dir,'norm');
        %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
        x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
        foundVerts(1,:) = [obj.G*x + obj.c]';

        % Find second (opposite) vertex
        dir = [0 -1 0];
        searchedDirs(2,:) = normalize(dir,'norm');
        %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
        x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
        foundVerts(2,:) = [obj.G*x + obj.c]';
    end

    if (foundVerts(1,:)-foundVerts(2,:) <= 1e-12) % Same vertex
        % Find first vertex
        dir = [0 0 1];
        searchedDirs(1,:) = normalize(dir,'norm');
        %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
        x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
        foundVerts(1,:) = [obj.G*x + obj.c]';

        % Find second (opposite) vertex
        dir = [0 0 -1];
        searchedDirs(2,:) = normalize(dir,'norm');
        %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
        x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
        foundVerts(2,:) = [obj.G*x + obj.c]';
    end


    if (foundVerts(1,:)-foundVerts(2,:) <= 1e-12) % Single point
        v = foundVerts(1,:);
        f = 1;
        return
    end

    % Continue if set is not a single point
    % Find 3rd vertex
    nullDirs = null(foundVerts(1,:)-foundVerts(2,:))';
    dir = nullDirs(1,:);
    %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
    x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
    extreme(1,:) = [obj.G*x + obj.c]';
    isNewVert(1) = sum(abs(nullDirs*(extreme(1,:)-foundVerts(1,:))'))>=1e-6;           % Tolerance
    dir = -dir;
    %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
    x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
    extreme(2,:) = [obj.G*x + obj.c]';
    isNewVert(2) = sum(abs(nullDirs*(extreme(2,:)-foundVerts(1,:))'))>=1e-6;           % Tolerance
    dir = nullDirs(2,:);
    %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
    x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
    extreme(3,:) = [obj.G*x + obj.c]';
    isNewVert(3) = sum(abs(nullDirs*(extreme(3,:)-foundVerts(1,:))'))>=1e-6;           % Tolerance
    dir = -dir;
    %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
    x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
    extreme(4,:) = [obj.G*x + obj.c]';
    isNewVert(4) = sum(abs(nullDirs*(extreme(4,:)-foundVerts(1,:))'))>=1e-6;           % Tolerance
    if max(isNewVert) == 0  % A line segment in 3D
        v = foundVerts;
        f = [1 2];
        return
    end

    % Continue if set is not a line segment in 3D
    indxPairs = [1 2; 1 3; 1 4; 2 3; 2 4; 3 4];
    coPlanar = zeros(6,1);
    for i = 1:6
        coPlanar(i) = abs(det([foundVerts(2,:)-foundVerts(1,:);extreme(indxPairs(i,1),:)-foundVerts(1,:);extreme(indxPairs(i,2),:)-foundVerts(1,:)])) <= 1e-6;
    end
    if min(coPlanar) == 1 % A planar set in 3D
        indx = find(isNewVert,1);
    %     basis = orth([foundVerts(2,:)-foundVerts(1,:);extreme(indx,:)-foundVerts(1,:)]')';
    %     reducedG = basis*obj.G;
        normVec = cross(foundVerts(2,:)-foundVerts(1,:),extreme(indx,:)-foundVerts(1,:));
        reducedG = obj.G(1:2,:); % Will not work if set is 'vertical';
        reducedObj = conZono(reducedG,obj.c(1:2),obj.A,obj.b);
        optPlot = plotOptions('Display','off');
        optPlot.SolverOpts = optSolver;
        [reducedV,~] = plot(reducedObj,optPlot);
    %     v = obj.c' + reducedV*basis;
        v = [reducedV 1/normVec(3)*(normVec*foundVerts(1,:)'-sum(normVec(1:2).*reducedV,2))];
        f = [1:size(v,1)];
        return
    end

    % Continue if set is not a 2D planar set in 3D
    % Form 3-simplex
    % indx = find(isNewVert,1);
    % foundVerts(3,:) = extreme(indx,:);
    % indx = find(coPlanar==0,1); 
    % foundVerts(4,:) = extreme(indx+1,:);
    indx = find(coPlanar==0,1);
    foundVerts(3,:) = extreme(indxPairs(indx,1),:);
    foundVerts(4,:) = extreme(indxPairs(indx,2),:);
    % Interior point
    interiorPoint = sum(foundVerts(1:4,:))/4;
    % Faces, vertices, and normals
    [f] = nchoosek(1:4,3);
    fNorms = zeros(0,3);
    for i = 1:size(f,1)
        [fN,fC] = faceNormal(foundVerts(f(i,1:3),:),interiorPoint);
        fNorms(i,:)= fN;
    end
    % Expand polytope
    maxIter = 1e6;                                      % Maximium iterations
    fIndex = 1;
    nVerts = 4;
    for i = 1:maxIter
        [dir,~] = faceNormal(foundVerts(f(fIndex,1:3),:),interiorPoint);
        %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
        x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
        extreme = [obj.G*x + obj.c]';
        [foundVerts,isNew] = addIfNew(foundVerts,extreme,dir);
        if isNew
            nVerts = nVerts + 1;
            % Find all faces that see the new vertex
            facesThatSee = findFacesThatSeeNewVertex(fNorms,foundVerts(f(:,1),:),extreme);
            sharedEdges = findSharedEdges(find(facesThatSee),f);
            facesToRemove = [];
            for indx = find(facesThatSee)
                removeCurrentFace = 0;
                if ismember(nVerts,f(indx,:))
                    continue
                end
                nVertsRow = find(f(indx,:),1,'last');
                startVerts = f(indx,1:nVertsRow);
                endVerts = [f(indx,2:nVertsRow) f(indx,1)];
                allCombos = [startVerts' endVerts'];
                allCombos = setdiff(sort(allCombos,2),sort(sharedEdges,2),'rows'); 
                newFaces = [allCombos nVerts*ones(size(allCombos,1),1)];
                minIndices = [];
                if isempty(newFaces)
                    removeCurrentFace = 1;
                end
                for j = 1:size(newFaces,1)
                    [fN,~] = faceNormal(foundVerts(newFaces(j,1:3),:),interiorPoint);
                    [minVal,minIndx] = min(sum(abs(fNorms-fN),2));
                    isNewFace = minVal>=1e-12;                    % Tolerance (Changed, 1e-6 too small => errors)
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
    f(f==0) = nan;
catch E
    if strcmp(E.identifier, 'PlotError:VertexNotFound')
        if checkEmpty(conZono(obj.G,obj.c,Aeq,beq))
            warning('zonoLAB:EmptyZonotope','Constrained zonotope is empty and cannot be plotted.')
            v = []; f = [];
            return
        end
        rethrow(E);
    else
        rethrow(E);
    end
end

end

% Local Functions
function [x] = findVertex(dir,G,Aeq,beq,lb,ub,optSolver)
    [x,~,~] = solveLP(dir*G,[],[],Aeq,beq,lb,ub,optSolver);
    if isnan(x)
        error('PlotError:VertexNotFound','Could not find a solution for a vertex while plotting')
    end
end

function [foundVerts,isNew] = addIfNew(foundVerts,extreme,dir)
    isNew = 0;
    % Check if this is the first found vertex
    if isempty(foundVerts)
        foundVerts = [foundVerts;extreme];
        isNew = 1;
        return
    end
    % Check if this is actually more extreme in the dir direction
    maxInDir = max(dir*foundVerts');
    if dir*extreme' - maxInDir <= 1e-6 % Threshold
        return % Not more extreme than existing vertices
    end
    % Check if this is the same as a vertex already found
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
