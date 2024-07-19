%% Creating Sets
% Creating a zonotope:
G = [1 0 1;0 1 1];  % Generator matrix with 3 generators
c = zeros(2,1);     % Center at the origin
Z = zono(G,c)       % Creates a zonotope

% Creating a constrained zonotope:
G = [1 0 1;0 1 1];   % Generator matrix with 3 generators
c = zeros(2,1);      % Center at the origin
A = ones(1,3);       % Constraint matrix with 1 constraint
b = 1;               % Constraint offset vector
C = conZono(G,c,A,b) % Creates a constrained zonotope

% Creating a hybrid zonotope:
Gc = [1 0 1;0 1 1];          % Continuous generator matrix with 3 generators
Gb = [1 2;2 1];              % Binary generator matrix with 2 generators
c = zeros(2,1);              % Center at the origin
Ac = ones(1,3);              % Continuous constraint matrix with 1 constraint
Ab = ones(1,2);              % Binary constraint matrix with 1 constraint
b = 1;                       % Constraint offset vector
H = hybZono(Gc,Gb,c,Ac,Ab,b) % Creates a hybrid zonotope

%% Set Operations
% Linear mapping of a zonotope from 2D to 3D using overloaded * (mtimes) operator:
M = [1 0;0 1;1 1];  % 3x2 matrix
Zm = M*Z            % Linear mapping

% Minkowski sum of a zonotope and constrained zonotope using overloaded + (plus) operator:
ZC = Z + C   % Minkowski sum

% Intersection of a hybrid zonotope and a constrained zonotope using the overloaded & (and) operator:
HC = H & C   % Intersection

%% Plotting
% Plot a random zonotope in 3D as transparent blue and output vertex and face information
seed = 2;
n = 3;
nG = 3;
Zr = randomSet(seed,'zono',n,nG);
figure;
[v,f] = plot(Zr,'b',0.1)

% Plot a hybrid zonotope with custom options
optsPlot = plotOptions;
optsPlot.FaceColor = 'r';
optsPlot.FaceAlpha = 0.5;
optsPlot.LineWidth = 3;
figure; 
plot(H,optsPlot)