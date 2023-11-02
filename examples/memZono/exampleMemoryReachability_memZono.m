clear; %close all;

%% Simulation Settings
N = 3; %< number of timesteps to do reachability
dt = 0.25;

%% System Definition
% Continous Time
A_ct = [0,1;-2,-1];
B_ct = [0;1];
C_ct = eye(2);
D_ct = 0;
sys_ct = ss(A_ct,B_ct,C_ct,D_ct);

% Discrete Time
sys = c2d(sys_ct,dt);
[A,B,C,D] = ssdata(sys);

% % Controller
% K_eig = [0.5;0.4];
% K = -place(A,B,K_eig);
% 
% % Estimator
% L_eig = 0.8*[1,1];
% L = place(A',C',L_eig);

% % Noise (not used? idk... it's weird)
% eta_ = 1; % \eta \in [-eta_,eta_]
% Eta_ = memZono()

% Initial and Final Sets %%%% (really need to get constructor better)
X_0 = memZono(zono(diag([1,2]),zeros(2,1)),'x_0');
X_F = memZono(zono(0.25*eye(2),ones(2,1)),'x_f');

% Input Set
U_nom = zono(1.5,0);
U_0 = memZono(U_nom,'u_0');

%% Reachability Calculation
% Initialization
X_all = [X_0; U_0]; %<--- cart prod but extendable to many of them
X_{1} = A*X_0 + B*U_0;
X_{1}.dimKeys = 'x_1';
U_{1} = memZono(U_nom,'u_1');
X_all = [X_all; X_{1}; U_{1}];
for k = 2:N
    % TODO: Add Measurement/Noise?
    % Next Time-step
    X_{k} = A*X_{k-1} + B*U_{k-1};
    X_{k}.dimKeys = sprintf('x_%d',k);
    % Input Calculation
    U_{k} = memZono(U_nom,sprintf('u_%d',k));
    % Name and save/label
    X_all = [X_all; X_{k}; U_{k}]; %<--- could add measurements/noise/etc too
end

%% Intersection
X_F.dimKeys = sprintf('x_%d',k-1);
Xall_inter_F = X_all & X_F; % <--- intersect common dimensions

%% Plotting
fig = figure;

% State plots
subplot(1,2,1);
hold on;
plot(X_F, 'all', 'g', 1);
plot(X_0, 'all', selectColor(0), 0.6)
plot(Xall_inter_F, sprintf('x_%d',0), selectColor(0), 0.6);
drawnow;
for k = 1:N-1
    plot(Xall_inter_F, sprintf('x_%d',k), selectColor(k), 0.6);
    plot(X_{k}, 'all', selectColor(k), 0.2);
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
plot([U_0; U_{1}],'all','b',0.2);
plot(Xall_inter_F,{'u_0','u_1'},'b',0.6);
drawnow
hold off;

axis equal;
xlim([-2 2]);
ylim([-2 2]);
xlabel('$u(1)$','Interpreter','latex');
ylabel('$u(2)$','Interpreter','latex');

%% local functions
function color = selectColor(i)
    colors = {'k','b','r'};
    color = colors{mod(i,length(colors))+1};
end