%% 1D
rng(1)
n = 1;
nG = 10;
G = 2*rand(n,nG)-1;
c = zeros(n,1);
Z = zono(G,c);

figure;
opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',1);
plot(Z,opt);

%% 2D
rng(1)
n = 2;
nG = 10;
% G = (2*rand(n,1)-1)*(2*rand(1,nG)-1); % 1D
G = 2*rand(n,nG)-1; % 2D
c = zeros(n,1);
Z = zono(G,c);

figure;
opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',1);
[v,f] = plot(Z,opt);

%% 3D
rng(1)
n = 3;
nG = 10;
% G = (2*rand(n,1)-1)*(2*rand(1,nG)-1); % 1D
% G = (2*rand(n,2)-1)*(2*rand(2,nG)-1); % 2D
% G =  [1 -1; 1 0; 1 1]*[1 2 3 4;2 3 7 9];
% G = [1 -1; 1 0; 1 -1]*[1 2;-1 2];
G = 2*rand(n,nG)-1; % 3D
c = zeros(n,1);
Z = zono(G,c);

figure;
tStart = tic;
opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',1);
[v,f] = plot(Z,opt);
drawnow
toc(tStart)

figure;
tStart = tic;
Box = Polyhedron('lb',-ones(nG,1),'ub',ones(nG,1));
P = c + G*Box;
plot(P,'color','r','alpha',1)
drawnow
toc(tStart)