clear; %close all;

%% Simulation Settings
N = 3; %< number of timesteps to do reachability
dt = 0.25;

%% System Definition
% DT-LTI System
A = [   1, 0.25;
     -0.5, 0.75];
B = [   0.025;
        0.25];

%% Reachability Setup
X_0 = zono(diag([1,2]),zeros(2,1));
X_F = zono(0.25*eye(2),ones(2.,1));
U_nom = zono(1.5,0);

% Dimensions
n = size(A,1);
p = size(B,2);
% q = size(C,1);

% Indexs
rx = {}; %<--- states
ru = {}; %<--- inputs
ri = 0; %<-- last index used

%% Reachability
% Initial Conditions
X_{1} = X_0; %<-- initial conditions
X_all = X_{1}; rx{1} = ri + (1:n); ri = ri + n;

% Time-evolution
for k = 1:N-1
    % Current Input
    U_{k} = U_nom;
    X_all = cartProd(X_all,U_{k}); ru{k} = ri + (1:p); ri = ri + p;

    % Step Update
    X_{k+1} = A*X_{k} + B*U_{k};
    X_all = extend_zonotope(X_all,X_{k+1}); rx{k+1} = ri + (1:n); ri = ri + n;
end

%% Intersection
R = zeros(n,X_all.n);
R(:,rx{N}) = eye(n);
X_inter = and(X_all,X_F,R);

%% Plotting
fig = figure;

% State plots
subplot(1,2,1);
hold on;
plot(X_F, 'g', 1);
drawnow;
for k = 1:N
    R = zeros(n,X_inter.n); R(:,rx{k}) = eye(n);
    plot(R*X_inter, selectColor(k), 0.6);
    plot(X_{k}, selectColor(k), 0.2);
    drawnow;
end
hold off;

axis equal;
xlim([-3 3]);
ylim([-3 3]);
xlabel('$x_1$','Interpreter','latex');
ylabel('$x_2$','Interpreter','latex');


% Input Plots
subplot(1,2,2);
hold on;
R = zeros(2,X_inter.n); R(1,ru{2}) = 1; R(2,ru{3}) = 1; %<--- figure out why this is weird
plot(R*X_inter,'b',0.6);
plot(cartProd(U_{2}, U_{3}),'b',0.2);
drawnow
hold off;

axis equal;
xlim([-2 2]);
ylim([-2 2]);
xlabel('$u(1)$','Interpreter','latex');
ylabel('$u(2)$','Interpreter','latex');

%% Local functions
function color = selectColor(i)
    colors = {'k','b','r'};
    color = colors{mod(i,length(colors))+1};
end

%% Functions

function [Z] = extend_zonotope(X_,Y)
% Extending X_ with Y. 
nx = size(X_.G,1);
ny = size(Y.G,1);
ncx = size(X_.A,1);
ncy = size(Y.A,1);
r_ngc_diff = 0;
l_ngc_diff = 0;
if size(Y.G,2) >= size(X_.G,2)
    %If X_ has f factors, then the first f factors of Y are the same as X_ 
    r_ngc_diff = size(Y.G,2) - size(X_.G,2); % num new cont. factors
elseif size(X_.G,2) > size(Y.G,2)
    %If Y has f factors, then the last f factors of X_ are the same as Y 
    l_ngc_diff = size(X_.G,2) - size(Y.G,2); % num new cont. factors
else
    error(sprintf('Ambiguous extension. One set should have more generators than the other (both binary and continuous).')) 
end

Z = conZono( [X_.G , zeros(nx,r_ngc_diff) ; zeros(ny,l_ngc_diff) , Y.G, ], ...
             [X_.c ; Y.c], ...
             [X_.A , zeros(ncx,r_ngc_diff) ; zeros(ncy,l_ngc_diff) , Y.A], ...
             [X_.b ; Y.b]);
end

% function [Z] = cartProd(X_,Y)
% % Stacking two independent zonotopes. 
% nx = size(X_.G,1);
% ny = size(Y.G,1);
% ncx = size(X_.A,1);
% ncy = size(Y.A,1);
% 
% Z = conZono( [X_.G , zeros(nx,size(Y.G,2)) ; zeros(ny,size(X_.G,2)) , Y.G, ], ...
%              [X_.c ; Y.c], ...
%              [X_.A , zeros(ncx,size(Y.G,2)) ; zeros(ncy,size(X_.G,2)) , Y.A], ...
%              [X_.b ; Y.b]);
% end