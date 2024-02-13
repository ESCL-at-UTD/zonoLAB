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

% Terminal Set
X_F = memZono(X_F,sprintf('x_%d',N));

switch 'affine&intersect' %'withOverload' 'affine&intersect'
    % case 'withOverload'
    %     % Time-evolution
    %     for k = 1:N-1
    %         % Current Input
    %         U_{k} = memZono(U_nom,sprintf('u_%d',k));

    %         % Step Update
    %         X_{k+1} = A*X_{k} + U_{k}.affine([],B,{},sprintf('x_%d',k));
    %         X_{k+1}.dimKeys = sprintf('x_%d',k+1);

    %         % Save Data
    %         % X_all = [X_all; U_{k}; X_{k+1}]; %<---- vertcat() = intersect()
    %         X_all = X_all & U_{k} & X_{k+1}; %<--- w/ and()

    %     end
    %     X_inter = X_all & X_F; % <--- intersect common dimensions
    case 'affine&intersect'
        % Time-evolution
        for k = 1:N-1
            % Current Input
            U_{k} = memZono(U_nom,sprintf('u_%d',k));

            % Step Update
            newDims = {sprintf('x_%d_1',k+1),sprintf('x_%d_2',k+1)};
            X_{k+1} = X_{k}.transform(U_{k}.transform([],B,{},newDims),A,{},newDims);

            % Save Data
            X_all = X_all.merge(U_{k});
            X_all = X_all.merge(X_{k+1});
        end
        X_inter = X_all.merge(X_F,'terminal_cons'); % <--- intersect common dimensions
end

%% Plotting
fig = figure;

% State plots
subplot(1,2,1);
hold on;
plot(X_F, 'all', 'g', 1);
drawnow;
for k = 1:N
    plot(X_inter, {sprintf('x_%d_1',k),sprintf('x_%d_2',k)}, selectColor(k), 0.6);
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