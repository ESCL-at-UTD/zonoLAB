% SLAM application of memZono class
clear; clc;

%% Plot settings 
set(0,'defaultLineLineWidth', 2)
set(0,'defaultAxesFontName' , 'Times')
set(0,'defaultTextFontName' , 'Times')
set(0,'defaultAxesFontSize' , 18)
set(0,'defaultTextFontSize' , 18)
set(0,'defaulttextinterpreter','latex')
set(0,'defaultlegendinterpreter','latex')
set(0,'defaultAxesGridLineStyle','-.')

%% System Definition

T = 12;         % number of time steps
num_l = 30;     % number of landmarks
theta = pi/6;   % rad between each measurement
F = [cos(theta),sin(theta);-sin(theta),cos(theta)]; % Rotation Matrix
G = eye(2);     % acceleration uncertainty multiplier
n = size(F,1);  % size
fig = 0;        % figure number

% initial conditions
X_nom = cell(1,T);              % nominal state zonotope (no constraints)
X = cell(1,T);                  % state zonotope with constraints
x0 = [5;0];                     % initial state
X0 = zono(0.01*eye(n),x0);      % initial state zonotope uncertainty
X_nom{1} = memZono(X0,'x_0e');  % state zonotope initial condition (labeled)
X_nom{1}.dimKeys = 'x_0e';      % label state zonotope initial dimensions (time-step zero)
X{1} = X_nom{1};                % state-landmark zonotope initial condition (labeled)
x = zeros(n,T);                 % initialize array for actual state values to exist in
x(:,1) = x0;                    % initial state
u = zeros(n,T);                 % inputs

V = zono(0.05*eye(2),zeros(2,1));   % system uncertainty
H = zono([0.075,0;0,0.05],[0;0]);   % measurement uncertainty

% landmark generation
landmarks = cell(1,num_l);      % exact landmarks
landmark_radius = 4;            % how far can the vehicle see
L = cell(T,num_l);              % all landmark location zonotopes
r_m = cell(T,num_l);            % all measurements of landmark distances
for i=1:num_l
    theta = 2*pi*rand();
    landmarks{i} = (2*rand()+2)*[cos(theta);sin(theta)];
end


%% simulation

for k = 1:T     % loop over every time step 

    fprintf('--- k = %i -----------------\n',k-1)   % print time step

    if k > 1                                                        % not initial condition
        % u(:,k-1) = random_sample_zonotope(U);                     % determine actual input
        mV = memZono(V,sprintf('v_%ie',k-1));                       % memZono version of the system uncertainty
        newDimKeys = memZono.genKeys(sprintf('x_%ie',k-1),1:2);                          % create dimension keys for the new time step
        X_nom{k} = X_nom{k-1}.transform([],F,{},newDimKeys) + mV.transform([],eye(2),{},newDimKeys);  % update state set factors
        %X_nom{k}.c = X_nom{k}.c + G*u(:,k);                         % update state set center
        x(:,k) = F*x(:,k-1) + G*u(:,k) + random_sample_zonotope(V); % determine actual state value
        X{k} = [X{k-1}; X_nom{k}];                                  % insert next time step to state-landmark zonotope
    end

    for i=1:num_l  % loop over each landmark

        r = landmarks{i}-x(:,k);    % actual vector displacement between vehicle and landmark
        theta = atan2(r(2),r(1));   % angle made relative to vehicle heading

        if norm(r) <= landmark_radius               % vehicle can see landmark

            fprintf('k = %i --> Landmark %i\n',k-1,i);  % print time-step and landmark measured
            RH = rotate_zonotope(H,theta);              % rotated measurement noise zonotope
            r_m{k,i} = r + random_sample_zonotope(RH);  % noisy measurement of landmark
            L{k,i} = RH + r_m{k,i};                     % set of landmark i from measurement at time-step k
            L{k,i} = memZono(L{k,i},['L_',num2str(k-1),'_',num2str(i)]);    % label factors appropriately
            L{k,i}.dimKeys = sprintf('x_%ie',k-1);                          % to ensure proper minkowski sum
            L{k,i} = plus(X_nom{k},L{k,i});                                 % add state uncertainty without constraints
            L{k,i}.dimKeys = sprintf('L_%ie',i);                            % label dimensions properly
            % add new landmark to zono by cartprod or generalized intersection
            X{k} = X{k}.merge(L{k,i},sprintf('L%i_k%i',i,k));
            % X{k} = [X{k}; L{k,i}];
        end
    end
end


%% plotting

for k = 1:T % loop over every time step

    % open figure
    if mod(k,2) ~= 0
        fig = fig + 1;  % new figure number
        figure(fig)     % open plot
        subplot(1,2,1)  % left subplot
        cla             % clear current figure
        grid on         % use grid
        box on          % create outline of plot
        axis equal      % scaling of axes is same
        hold on         % allow multiple plots
    else
        subplot(1,2,2)  % right subplot
        cla             % clear current figure
        grid on         % use grid
        box on          % create outline of plot
        axis equal      % scaling of axes is same
        hold on         % allow multiple plots
    end
    title(sprintf('$k = %i$',k-1),'interpreter','Latex');

    % plot landmarks
    for i = 1:num_l

        % find the last time-step when this landmark was measured
        last_measurement = k;   % set time of last measurement to the current time-step
        while isempty(L{last_measurement,i})    % only proceed when a measurement exists
            % go to previous time-step to check for measurements
            last_measurement = last_measurement - 1;
            % break if no measurement exists
            if last_measurement < 1
                break
            end
        end

        % plot measurements of current time-step
        if last_measurement == k
            plot(r_m{k,i}(1)+x(1,k),r_m{k,i}(2)+x(2,k),'b.','MarkerSize',8)
            plot([x(1,k),r_m{k,i}(1)+x(1,k)],[x(2,k),r_m{k,i}(2)+x(2,k)],'b','LineWidth',0.5);
            drawnow;
        end

        % plot zonotopes
        if last_measurement >= 1                        % ensure landmark has been measured
            plot(X{k},X{k}.keysStartsWith(sprintf('L_%ie',i)).dimKeys,'r',0.6);    % plot landmark zonotope
            plot(landmarks{i}(1),landmarks{i}(2),'bx'); % plot landmarks accurate location
            drawnow;
        end       
    end

    % plot state zonotopes
    for prev_k = 1:k                                        % plot all previous time step states
        plot(x(1,prev_k),x(2,prev_k),'k.','MarkerSize',12); % plot actual state
        plot(X_nom{prev_k},'all','k',0.2);                  % plot unconstrained state zonotope
        plot(X{k},X{k}.keysStartsWith(sprintf('x_%ie',prev_k-1)).dimKeys,'g',0.3);     % plot constrained state zonotope
        drawnow;
    end
    % xlim([-6 6]);ylim([-6 6]);
    xlim([-5.75 5.75]);ylim([-5.75 5.75]);
end

% exportgraphics(gcf,'SLAM_ZONO.pdf','ContentType','vector')


%%

function Z_prime = rotate_zonotope(Z, radians)
    center = Z.c;
    Z.c = Z.c - center;
    R = [cos(radians) -sin(radians); sin(radians) cos(radians)];
    Z_prime = R * Z;
    Z_prime.c = Z_prime.c + center;
end

% random_sample_zonotope()
function s = random_sample_zonotope(z)
    g = 2*rand([z.nG,1])-1;
    s = z.c + z.G*g;
end
