clear;

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

%% Reachability Calculation
% Initial Conditions
X_{1} = memZono(X_0,'x_1');
X_all = X_{1};

% Time-evolution
for k = 1:N
    % Current Input
    U_{k} = memZono(U_nom,sprintf('u_%d',k));

    % Step Update
    X_{k+1} = A*X_{k} + B*U_{k};
    X_{k+1}.dimKeys = sprintf('x_%d',k+1);

    % Save Data
    X_all = [X_all; U_{k}; X_{k+1}];
end

%% Intersection
X_F = memZono(X_F,sprintf('x_%d',k+1));
X_inter = X_all & X_F; % <--- intersect common dimensions

%% Plotting
fig = figure;

% State plots
subplot(1,2,1);
hold on;
plot(X_F, 'all', 'g', 1);
drawnow;
for k = 1:N
    plot(X_inter, sprintf('x_%d',k), selectColor(k), 0.6);
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
plot([U_{1}; U_{2}],'all','b',0.2);
plot(X_inter,{'u_1','u_2'},'b',0.6);
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