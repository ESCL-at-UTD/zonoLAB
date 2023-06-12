function [v,f] = plotConZono2D(obj,opts)

% Standardized header

Aeq = sparse(obj.A);
beq = [obj.b];
lb = -ones(obj.nG,1);
ub =  ones(obj.nG,1);

% Find first vertex
dir = [1 0];
searchedDirs(1,:) = normalize(dir,'norm');
[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,opts);
foundVerts(1,:) = [obj.G*x + obj.c]';

% Find second (opposite) vertex
dir = [-1 0];
searchedDirs(2,:) = normalize(dir,'norm');
[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,opts);
foundVerts(2,:) = [obj.G*x + obj.c]';

if (foundVerts(1,:)-foundVerts(2,:) <= 1e-12) % Same vertex
    % Find first vertex
    dir = [0 1];
    searchedDirs(1,:) = normalize(dir,'norm');
    [x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,opts);
    foundVerts(1,:) = [obj.G*x + obj.c]';

    % Find second (opposite) vertex
    dir = [0 -1];
    searchedDirs(2,:) = normalize(dir,'norm');
    [x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,opts);
    foundVerts(2,:) = [obj.G*x + obj.c]';
end

if (foundVerts(1,:)-foundVerts(2,:) <= 1e-12) % Single point
    v = foundVerts(1,:);
    f = 1;
    return
end

vertDiff = foundVerts(1,:) - foundVerts(2,:);
dir = [-vertDiff(2) vertDiff(1)];
[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,opts);
extreme = [obj.G*x + obj.c]';
isNewVert(1) = dir*extreme' >= 1e-6;                        % Tolerance
dir = -dir;
[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,opts);
extreme = [obj.G*x + obj.c]';
isNewVert(2) = dir*extreme' >= 1e-6;                        % Tolerance
if max(isNewVert) == 0  % A line segment in 2D
    v = foundVerts;
    f = [1 2];
    return
end

firstVert = foundVerts(1,:);
endVert = foundVerts(2,:);
% Compute center
centerVert = (firstVert + endVert)/2;
% Compute vertex angles
vertDiff = firstVert - centerVert;
referenceAngle = atan2(vertDiff(:,2),vertDiff(:,1));
foundVertDirs(1) = 0;
% Farthest CW vert
foundVertDirs(2) = -pi;
treeStruct(1,3) = 2;
treeStruct(2,:) = [1 zeros(1,2)];
% Farthest CW vert
foundVerts(3,:) = endVert;
foundVertDirs(3) = pi;
treeStruct(1,2) = 3;
treeStruct(3,:) = [1 zeros(1,2)];
% Initialize number of vertices
nVerts = 3;
% Order the vertices
[vertOrder,~] = listInOrder(treeStruct,1,zeros(nVerts,1),1);
% Find the remaining vertices
maxIter = 1e6;
indx = 1;
for i = 1:maxIter
    vertA = foundVerts(vertOrder(indx),:);
    vertB = foundVerts(vertOrder(indx+1),:);
    vertDiff = vertA - vertB;
    dir = -[-vertDiff(2) vertDiff(1)];
    isNewdir = min(sum(abs(normalize(dir,'norm')-searchedDirs),2))>=1e-6;   % Tolerance
    if isNewdir
        searchedDirs = [searchedDirs;normalize(dir,'norm')];                % Growing list
        [x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,opts);
        extreme = [obj.G*x + obj.c]';
        vertDiffs = extreme - centerVert;
        vertAngle = atan2(vertDiffs(:,2),vertDiffs(:,1));
        newVertDir = vertAngle-referenceAngle;
        newVertDir(newVertDir<-pi) = newVertDir(newVertDir<-pi) + 2*pi;
        newVertDir(newVertDir>pi) = newVertDir(newVertDir>pi) - 2*pi;
        orderDirs = foundVertDirs(vertOrder);
        isNew = min(abs(newVertDir-orderDirs))>=1e-6;                   % Tolerance
    else
        isNew = 0;
    end
    if isNew
        nVerts = nVerts + 1;
        foundVerts(nVerts,:) = extreme;
        foundVertDirs(nVerts) = newVertDir;
        if treeStruct(vertOrder(indx),3) == 0
            treeStruct(vertOrder(indx),3) = nVerts;
            treeStruct(nVerts,:) = [vertOrder(indx) 0 0];
        else
            treeStruct(vertOrder(indx+1),2) = nVerts;
            treeStruct(nVerts,:) = [vertOrder(indx+1) 0 0];
        end
        [vertOrder,~] = listInOrder(treeStruct,1,zeros(nVerts,1),1);
    else
        if vertOrder(indx+1) == 2
            break
        else
            indx = indx + 1;
        end
    end
end

nVerts = nVerts - 1; % Added to remove extra vertex
v = foundVerts(vertOrder,:);
v = v(1:end-1,:); % Added to remove extra vertex
f = [1:nVerts];

end

% Local functions
function [vertOrder,indx] = listInOrder(treeStruct,row,vertOrder,indx)
    if treeStruct(row,2) ~= 0
        [vertOrder,indx] = listInOrder(treeStruct,treeStruct(row,2),vertOrder,indx);
    end
    vertOrder(indx) = row;
    indx = indx + 1;
    if treeStruct(row,3) ~= 0
        [vertOrder,indx] = listInOrder(treeStruct,treeStruct(row,3),vertOrder,indx);
    end   
end