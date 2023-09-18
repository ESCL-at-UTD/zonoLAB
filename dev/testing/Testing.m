% c = zeros(2,1);
% G = eye(2);
% Z = zono(G,c);
% Y = Z;
% ZAll(1) = Z;
% ZAll(2) = Y;
% ZAll
% 
% tic
% X = Z + Y;
% toc
% % X
% 
% A = ones(1,2);
% b = 1;
% C = conZono(G,c,A,b);
% CAll(1) = C;
% CAll(2) = C;
% CAll
% H = hybZono(G,G,c,A,A,b);
% 
% HAll(1) = H;
% HAll(2) = H;
% HAll
%%
rng(2)
nG = 10;
G = 2*rand(2,nG)-1;
c = zeros(2,1);
Z4 = zono(G,c);

figure
opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',0.);
plot(Z4,opt);
% plot(Z4,'b',0.1)

% CORA_Zono = zonotope(c,G);
% figure;
% tic;
% plot(CORA_Zono)
% toc


% figure;
% tic;
% plot(Z4)
% toc

%%

% figure;
% tic;
% Box = Polyhedron('lb',-ones(Z4.nG,1),'ub',ones(Z4.nG,1));
% P = plus(Z4.c,affineMap(Box,Z4.G));
% plot(P,'color','r','alpha',1)
% toc
%%
rng(2)
nG = 400;
G = 2*rand(3,nG)-1;
c = zeros(3,1);
Z4 = zono(G,c);
tStart = tic;
figure;
plot(Z4,'r',0.1)
% twoCombos = ff2n(2);
% twoCombos(twoCombos == 0) = -1;
% twoCombos(4:-1:3,:) = twoCombos(3:4,:);
% allCombos = nchoosek(1:nG,2);
% V = zeros(8*size(allCombos,1),3);
% % figure; hold on
% for i = 1:size(allCombos,1)
%     verts = twoCombos*G(:,allCombos(i,:))';
%     FN = cross(G(:,allCombos(i,1)),G(:,allCombos(i,2)));
%     remainingGens = setdiff([1:nG],allCombos(i,:));
%     dotProds = FN'*G(:,remainingGens);
%     signs = sign(dotProds);
%     centerA = signs*G(:,remainingGens)';
%     centerB = -signs*G(:,remainingGens)';
%     V([1:4]+(i-1)*8,:) = centerA + verts;
%     V([5:8]+(i-1)*8,:) = centerB + verts;
% end
% F = reshape([1:size(V,1)],4,[])';
% figure;hold on
% patch('Faces',F,'Vertices',V,'FaceColor','red','FaceAlpha',0.1)
% view(3)
% drawnow
toc(tStart)

% figure;
% tStart = tic;
% Box = Polyhedron('lb',-ones(Z4.nG,1),'ub',ones(Z4.nG,1));
% P = plus(Z4.c,affineMap(Box,Z4.G));
% plot(P,'color','r','alpha',1)
% drawnow
% toc(tStart)
%
% This works but is very inefficient
% allCombos = ff2n(nG);
% allCombos(allCombos == 0) = -1;
% V = allCombos*G';
% F = convhull(V);
% TR = triangulation(F,V);
% centers = incenter(TR);
% FN = faceNormal(TR);
% figure;hold on
% patch('Faces',F,'Vertices',V,'FaceColor','red','FaceAlpha',0.1)
% quiver3(centers(:,1),centers(:,2),centers(:,3), ...
%      FN(:,1),FN(:,2),FN(:,3),0.5,'color','r');
% view(3)
% drawnow
% This is intended to combine triangulations
% F_new = zeros(size(FN,1)/2,4);
% indxList = 1:size(FN,1);
% alignmentMatrix = FN*FN';
% [row column] = find(abs(alignmentMatrix - 1) <= 1e-12);
% M = [row column];
% M(M(:,1)==M(:,2),:) = [];
% M(M(:,1)>=M(:,2),:) = [];
% F_new =  [F(M(:,1),:) F(M(:,2),:)];
% sortedF_new = sort(F_new,2);
% [row column] = find([diff(sortedF_new,1,2)~=0]');
% ind = sub2ind(size(sortedF_new'),row+1,column);
% F_new2 = [F_new(:,1) reshape(sortedF_new_T(ind),3,[])'];
% figure;hold on
% patch('Faces',F_new3,'Vertices',V,'FaceColor','red','FaceAlpha',0.1)
% view(3)
% drawnow
% toc(tStart)
% figure
% opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',0.);
% plot(Z4,opt);


%%
% Define set 1 dimensions and matrices
n1 = 2;
nGc1 = 5;
nGb1 = 3;
nC1 = 2;
Z1  = randomSet(1,'zono',n1,nGc1,nGb1,nC1);
Zc1 = randomSet(1,'conZono',n1,nGc1,nGb1,nC1);
Zh1 = randomSet(1,'hybZono',n1,nGc1,nGb1,nC1);
% Define set 2 dimensions and matrices
n2 = 3;
nGc2 = 8;
nGb2 = 3;
nC2 = 3;
Z2  = randomSet(2,'zono',n2,nGc2,nGb2,nC2);
Zc2 = randomSet(2,'conZono',n2,nGc2,nGb2,nC2);
Zh2 = randomSet(2,'hybZono',n2,nGc2,nGb2,nC2);

rng(1)
R = rand(n2,n1);
figure; hold on
plot(Z1,'r',0.1)
% plot(Z2,'b',0.1)
plot(and(Z1,Z2,R),'g',0.5)

figure; hold on
plot(Zc1,'r',0.1)
% plot(Zc2,'b',0.1)
plot(and(Zc1,Zc2,R),'g',0.5)

figure; hold on
plot(Zh1,'r',0.1)
% plot(Zh2,'b',0.1)
plot(and(Zh1,Zh2,R),'g',0.5)

% figure; hold on
% plot(Z1,'r',0.1)
% plot(Z2,'b',0.1)
% plot(Z1&Z2,'g',0.5)
% 
% figure; hold on
% plot(Zc1,'r',0.1)
% plot(Zc2,'b',0.1)
% plot(Zc1&Zc2,'g',0.5)
% 
% figure; hold on
% plot(Zh1,'r',0.1)
% plot(Zh2,'b',0.1)
% plot(Zh1&Zh2,'g',0.5)

%% Convex Hull
n = 2;
nG = 5;
nC = ceil(nG/2);
X  = randomSet(1,'conZono',n,nG,[],nC);
Y  = randomSet(2,'conZono',n,nG,[],nC);
Z = convexHull(X,Y);
figure; hold on
tStart = tic;
[v,f] = plot(X,'b',0.1);
[v,f] = plot(Y,'m',0.1);
[v,f] = plot(Z,'g',0.1);
toc(tStart)


%% Pontryagin Difference
n = 2;
nG = 7;
nC = ceil(nG/2);
nb = 2;
X  = 30*randomSet(1,'hybZono',n,nG,nb,nC);
Y  = randomSet(2,'zono',n,nG,[],[]);
Z = pontryDiff(X,Y);
figure; hold on
tStart = tic;
[v,f] = plot(X,'b',0.1);
[v,f] = plot(Y,'m',0.1);
[v,f] = plot(Z,'g',0.1);
toc(tStart)