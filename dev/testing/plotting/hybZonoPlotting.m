%% 1D
rng(1)
n = 1;
nGc = 10;
nGb = 2;
nC = 1;
Gc = 2*rand(n,nGc)-1;
Gb = 20*rand(n,nGb)-1;
c = zeros(n,1);
Ac = 2*rand(nC,nGc)-1;
Ab = 2*rand(nC,nGb)-1;
b = 2*rand(nC,1)-1;
% Ac = zeros(0,nGc);
% Ab = zeros(0,nGb);
% b = [];
Z = hybZono(Gc,Gb,c,Ac,Ab,b);

optSolver = solverOptions('milpSolver','gurobi');

figure;
tStart = tic;
opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',1,'SolverOpts',optSolver);
[v,f] = plot(Z,opt);
drawnow
toc(tStart)

figure;hold on
tStart = tic;
[leaves] = getLeaves(Z,optSolver);
for i = 1:size(leaves,2)
    Box = Polyhedron('lb',-ones(nGc,1),'ub',ones(nGc,1),'He',[Ac b-Ab*leaves(:,i)]);
    P = (c + Gb*leaves(:,i)) + Gc*Box;
    plot(P,'color','r','alpha',1)
%     drawnow
end
toc(tStart)

%% 2D
rng(1)
n = 2;
nGc = 10;
nGb = 4;
nC = 1;
% Gc = (2*rand(n,1)-1)*(2*rand(1,nGc)-1); % 1D
Gc = 2*rand(n,nGc)-1;
Gb = 20*rand(n,nGb)-1;
c = zeros(n,1);
Ac = 2*rand(nC,nGc)-1;
Ab = 2*rand(nC,nGb)-1;
b = 2*rand(nC,1)-1;
Z = hybZono(Gc,Gb,c,Ac,Ab,b);

optSolver = solverOptions('milpSolver','gurobi');

figure;
tStart = tic;
opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',1,'SolverOpts',optSolver);
[v,f] = plot(Z,opt);
drawnow
toc(tStart)

figure;hold on
tStart = tic;
[leaves] = getLeaves(Z,optSolver);
for i = 1:size(leaves,2)
    Box = Polyhedron('lb',-ones(nGc,1),'ub',ones(nGc,1),'He',[Ac b-Ab*leaves(:,i)]);
    P = (c + Gb*leaves(:,i)) + Gc*Box;
    plot(P,'color','r','alpha',1)
%     drawnow
end
toc(tStart)

%% 3D
% profile on
rng(1)
n = 3;
nGc = 10;
nGb = 2;
nC = 1;
% Gc = (2*rand(n,1)-1)*(2*rand(1,nGc)-1); % 1D
% Gc = (2*rand(n,2)-1)*(2*rand(2,nGc)-1); % 2D
Gc = 2*rand(n,nGc)-1;
Gb = 20*rand(n,nGb)-1;
c = zeros(n,1);
Ac = 2*rand(nC,nGc)-1;
Ab = 2*rand(nC,nGb)-1;
b = 2*rand(nC,1)-1;
Z = hybZono(Gc,Gb,c,Ac,Ab,b);

optSolver = solverOptions('milpSolver','gurobi');
optSolver.lpSolver = 'gurobi';
figure;
tStart = tic;
opt = plotOptions('FaceColor','b','FaceAlpha',0.1,'EdgeAlpha',1,'SolverOpts',optSolver);
[v,f] = plot(Z,opt);
drawnow
toc(tStart)

figure;hold on
tStart = tic;
[leaves] = getLeaves(Z,optSolver);
for i = 1:size(leaves,2)
    Box = Polyhedron('lb',-ones(nGc,1),'ub',ones(nGc,1),'He',[Ac b-Ab*leaves(:,i)]);
    P = (c + Gb*leaves(:,i)) + Gc*Box;
    plot(P,'color','r','alpha',1)
%     drawnow
end
toc(tStart)

% profile viewer