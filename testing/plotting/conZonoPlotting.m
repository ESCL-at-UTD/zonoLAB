%% 1D
rng(1)
n = 1;
nG = 10;
nC = 4;
G = 2*rand(n,nG)-1;
c = zeros(n,1);
A = 2*rand(nC,nG)-1;
b = 2*rand(nC,1)-1;
Z = conZono(G,c,A,b);

figure;
tStart = tic;
optSolver = solverOptions('lpSolver','linprog');
opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',1,'SolverOpts',optSolver);
[v,f] = plot(Z,opt);
drawnow
toc(tStart)

figure;
tStart = tic;
Box = Polyhedron('lb',-ones(nG,1),'ub',ones(nG,1),'He',[A b]);
P = c + G*Box;
plot(P,'color','r','alpha',1)
drawnow
toc(tStart)

%% 2D
rng(1)
n = 2;
nG = 10;
nC = 4;
G = (2*rand(n,1)-1)*(2*rand(1,nG)-1); % 1D
% G = 2*rand(n,nG)-1;
c = zeros(n,1);
A = 2*rand(nC,nG)-1;
b = 2*rand(nC,1)-1;
Z = conZono(G,c,A,b);

figure;
tStart = tic;
optSolver = solverOptions('lpSolver','gurobi');
opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',1,'SolverOpts',optSolver);
[v,f] = plot(Z,opt);
drawnow
toc(tStart)

figure;
tStart = tic;
Box = Polyhedron('lb',-ones(nG,1),'ub',ones(nG,1),'He',[A b]);
P = c + G*Box;
plot(P,'color','r','alpha',1)
drawnow
toc(tStart)

%% 3D
rng(1)
n = 3;
nG = 10;
nC = 4;
% G = (2*rand(n,1)-1)*(2*rand(1,nG)-1); % 1D
G = (2*rand(n,2)-1)*(2*rand(2,nG)-1); % 2D
% G = 2*rand(n,nG)-1;
c = zeros(n,1);
A = 2*rand(nC,nG)-1;
b = 2*rand(nC,1)-1;
Z = conZono(G,c,A,b);

figure;
tStart = tic;
optSolver = solverOptions('lpSolver','gurobi');
opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',1,'SolverOpts',optSolver);
% opt.Display = 'individual';
[v,f] = plot(Z,opt);
drawnow
toc(tStart)

figure;
tStart = tic;
Box = Polyhedron('lb',-ones(nG,1),'ub',ones(nG,1),'He',[A b]);
P = c + G*Box;
plot(P,'color','r','alpha',1)
drawnow
toc(tStart)
