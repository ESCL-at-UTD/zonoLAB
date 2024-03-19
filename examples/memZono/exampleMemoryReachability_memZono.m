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

n = size(A,1); m = size(B,2);

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

switch 'withOverload' %'withOverload' 'transform&merge' 'fullStack'
    case 'withOverload'
        % Time-evolution
        for k = 1:N-1
            % Current Input
            U_{k} = memZono(U_nom,sprintf('u_%d',k));

            % Step Update
            newDims = sprintf('x_%d',k+1);
            X_{k+1} = X_{k}.linMap(A,newDims) + U_{k}.linMap(B,newDims);

            % Save Data
            X_all = [
                X_all; 
                U_{k}; 
                X_{k+1}]; %<---- vertcat() = merge() [really just cartprod()]

        end
        X_inter = X_all.merge(X_F,'terminal_cons'); %<-- merge does the intersection

    case 'transform&merge'
        % Time-evolution
        for k = 1:N-1
            % Current Input
            U_{k} = memZono(U_nom,sprintf('u_%d',k));

            % Step Update
            newDims = {sprintf('x_%d_1',k+1),sprintf('x_%d_2',k+1)};
            X_{k+1} = X_{k}.transform(U_{k}.transform([],B,{},newDims),A,{},newDims); %<== transform has affine A x + B

            % Save Data
            X_all = X_all.merge(U_{k});
            X_all = X_all.merge(X_{k+1});
        end
        X_inter = X_all.merge(X_F,'terminal_cons'); % <--- intersect common dimensions


    case 'fullStack'
        % Label Functions
        xLabels = @(k) memZono.genKeys(sprintf('x_%d',k),1:n);
        % uLabels = @(k) memZono.genKeys(sprintf('u_%d',k),1:m);
        uLabels = @(k) {sprintf('u_%d',k)}; %<== 1D u_k

        % Time-Evolution
        for k = 1:N-1
            % Current Input
            X_all = X_all.merge(memZono(U_nom,uLabels(k)));

            % Step Update
            % switch 'indirectTransform' % 'directTransform' 'indirectTransform' 'combineTransform' 'overloading'
            %     case 'directTransform'
            % X_all = X_all.combine(...
            %     X_all.transform( ...
            %         transform(X_all(uLabels(k)),[],B,uLabels(k),xLabels(k)),...
            %         A, xLabels(k), xLabels(k+1))...
            %             );
            % 
            %     case 'combineTransform'
            % X_all = horzcat(X_all,...
            %     X_all.transform([],A,xLabels(k),xLabels(k+1)),...
            %         X_all.transform([],B,uLabels(k),xLabels(k+1))...
            % );
                % case 'overloading'
            
            % subsref(), linMap(), plus()
            X_all = [X_all;
                linMap(X_all(xLabels(k)),A,xLabels(k+1))... 
                    + linMap(X_all(uLabels(k)),B,xLabels(k+1));
            ];

            % end
        end

        X_inter = X_all.merge(X_F,'terminal_cons');

        for k = 1:N-1
            X_{k} = X_all(xLabels(k));
            U_{k} = X_all(uLabels(k));
        end
        X_{N} = X_all(xLabels(N));

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