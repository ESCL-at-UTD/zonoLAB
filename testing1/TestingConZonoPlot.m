rng(1)
nG = 40;
nC = 10;
G = 2*rand(2,nG)-1;
c = zeros(2,1);
A = 2*rand(nC,nG)-1;
b = 2*rand(nC,1)-1;
obj = conZono(G,c,A,b);

%% CORA plot
% CORA = conZonotope(c,G,A,b);
% figure;
% tStart = tic;
% plot(CORA)
% drawnow
% toc(tStart)
% 
% % MPT plot
% 
% figure;
% tStart = tic;
% Box = Polyhedron('lb',-ones(obj.nG,1),'ub',ones(obj.nG,1),'He',[obj.A obj.b]);
% P = obj.c + obj.G*Box;
% plot(P,'color','r','alpha',1)
% drawnow
% toc(tStart)

%%
figure; hold on
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
foundVerts = zeros(maxVerts*nLeaves,2);
vertOrder = zeros(maxVerts*nLeaves,1);
treeStruct = zeros(maxVerts*nLeaves,3);
searchedDirs = zeros(0,2);
totalVerts = 0;
foundVertDirs = zeros(maxVerts*nLeaves,1);
firstIndices = [1:maxVerts:maxVerts*nLeaves]';
% Find first vertex
dir = [1 0];
searchedDirs = [searchedDirs;normalize(dir,'norm')];
[extreme] = getVertices(obj,model,params,dir);
firstVert = extreme;
foundVerts(firstIndices,:) = firstVert;
% plot(extreme(:,1),extreme(:,2),'ok')
% Find opposite vertex
dir = -[1 0];
searchedDirs = [searchedDirs;normalize(dir,'norm')];
[extreme] = getVertices(obj,model,params,dir);
endVert = extreme;
% plot(extreme(:,1),extreme(:,2),'ok')
% Compute center
centerVert = (foundVerts(firstIndices,:) + endVert)/2;
% plot(centerVert(:,1),centerVert(:,2),'og')
% Compute vertex angles
vertDiffs = firstVert - centerVert;
referenceAngles = atan2(vertDiffs(:,2),vertDiffs(:,1));
foundVertDirs(firstIndices) = zeros(nLeaves,1);
% Farthest CW vert
foundVerts(firstIndices+1,:) = endVert;
foundVertDirs(firstIndices+1) = -pi;
treeStruct(firstIndices,3) = firstIndices+1;
treeStruct(firstIndices+1,:) = [firstIndices zeros(nLeaves,2)];
% Farthest CW vert
foundVerts(firstIndices+2,:) = endVert;
foundVertDirs(firstIndices+2) = pi;
treeStruct(firstIndices,2) = firstIndices+2;
treeStruct(firstIndices+2,:) = [firstIndices zeros(nLeaves,2)];
% Find the remaining vertices
nVerts = 3*ones(nLeaves,1);
for i = 1:nLeaves
    [vertOrder(firstIndices(i)+[0:nVerts(i)-1]),~] = listInOrder(treeStruct,firstIndices(i),zeros(nVerts(i),1),1);
end
maxIter = 1e6;
indx = 1;
j = 1;
for i = 1:maxIter
    vertA = foundVerts(vertOrder(indx),:);
    vertB = foundVerts(vertOrder(indx+1),:);
    vertDiff = vertA - vertB;
    dir = -[-vertDiff(2) vertDiff(1)];
    isNewdir = min(sum(abs(normalize(dir,'norm')-searchedDirs),2))>=1e-6;   % Tolerance
    if isNewdir
        searchedDirs = [searchedDirs;normalize(dir,'norm')];                % Growing list
        [extreme] = getVertices(obj,model,params,dir);
        vertDiffs = extreme - centerVert;
        vertAngles = atan2(vertDiffs(:,2),vertDiffs(:,1));
        newVertDirs = vertAngles-referenceAngles;
        newVertDirs(newVertDirs<-pi) = newVertDirs(newVertDirs<-pi) + 2*pi;
        newVertDirs(newVertDirs>pi) = newVertDirs(newVertDirs>pi) - 2*pi;
        orderDirs = foundVertDirs(vertOrder(firstIndices(j)+[0:nVerts(j)-1]));
        isNew = min(abs(newVertDirs(j)-orderDirs))>=1e-6;                   % Tolerance
    else
        isNew = 0;
    end

    if isNew
        nVerts(j) = nVerts(j) + 1;
        foundVerts(firstIndices(j)+nVerts(j)-1,:) = extreme(j,:);
        foundVertDirs(firstIndices(j)+nVerts(j)-1) = newVertDirs(j);
        if treeStruct(vertOrder(indx),3) == 0
            treeStruct(vertOrder(indx),3) = firstIndices(j)+nVerts(j)-1;
            treeStruct(firstIndices(j)+nVerts(j)-1,:) = [vertOrder(indx) 0 0];
        else
            treeStruct(vertOrder(indx+1),2) = firstIndices(j)+nVerts(j)-1;
            treeStruct(firstIndices(j)+nVerts(j)-1,:) = [vertOrder(indx+1) 0 0];
        end
        [vertOrder(firstIndices(j)+[0:nVerts(j)-1]),~] = listInOrder(treeStruct,firstIndices(j),zeros(nVerts(j),1),1);
%         plot(extreme(j,1),extreme(j,2),'ok')
    else
        if vertOrder(indx+1) == firstIndices(j)+1
            totalVerts(j) = sum(nVerts);
            if j == nLeaves
                break
            else
                j = j + 1;
                indx = firstIndices(j);
            end
        else
            indx = indx + 1;
        end
    end
end
nVerts = nVerts - 1; % Added to remove extra vertex (May not work with hybZono)
v = foundVerts(vertOrder(vertOrder~=0),:);
f = nan*ones(nLeaves,max(nVerts));
v = v(1:end-1,:); % Added to remove extra vertex (May not work with hybZono)
% f = f(:,1:end-1); % Added to remove extra vertex (May not work with hybZono)
count = 1;
for i2 = 1:nLeaves
    f(i2,1:nVerts(i2)) = count+[0:nVerts(i2)-1];
    count = count + nVerts(i2);
end
patch('Faces',f,'Vertices',v,'FaceColor','m')
drawnow
toc(tStart)

%%
function [extreme] = getVertices(obj,model,params,dir)
model.obj = dir*obj.G;
result = gurobi(model,params);
extreme = [obj.G*result.x + obj.c]';
end

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