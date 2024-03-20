clear; close all;
%% Trig Functions as Hybrid Zonotopes
x_max = pi/4;
x_min = -x_max;
n_points = 5; % Number of points (should be an odd number)

sinZono = makeSinX(n_points,[-pi/4 pi/4]);
cosZono = makeCosX(n_points,[-pi/4 pi/4]);

%% Using the Hybrid Zonotope approximation of Sin/Cos as the dynamics of a unicycle:
%     x_{k+1} = x_k + v*cos( theta_k )
%     y_{k+1} = y_k + v*sin( theta_k )
% theta_{k+1} = theta_k + u_k
% (the dt term could be absorbed into v and u_k.)

% initial state
Z0 = zono([0.01,0.02,0;0,0.01,0;0,0,5*pi/180],[0;0;0]);

N = 3; % number of time steps
v = 0.1; % constant forward velocity
U = zono([5/180*pi],5/180*pi); % angular velocity

% A naive non-memory approach has 9216 leaves at time step 3
% it also does not preserve the dependency across time, so it 
do_direct = 0;
if do_direct
    Z_direct = evolveUnicycle(sinZono,cosZono,v,U,Z0,N);
end

% A memory approach has 28 leaves at time step 3 and 208 at time step 4
Z = evolveUnicycleMemZono(sinZono,cosZono,v,U,Z0,N);

%% Plotting
% merge all time steps together
% to do so we need to relabel the dimKeys
% Z{i}({'x','y','theta'}) ensures that the dimKeys are in the expected order for relabeling
i = 1;
liftZ = copy(Z{i}({'x','y','theta'}),{sprintf('x_%i',i),sprintf('y_%i',i),sprintf('theta_%i',i)});
for i = 2:N
    liftZ = merge(liftZ, copy(Z{i}({'x','y','theta'}),{sprintf('x_%i',i),sprintf('y_%i',i),sprintf('theta_%i',i)}));
end

% imagine that new data is received that informs the y value at time step 3
Z_newdata = memZono(zono(0.005,0.025),'z_pin');
Z_newdata.dimKeys = 'y_3';
liftZ_newdata = merge(liftZ,Z_newdata,'pin_merge');

for i = 1:N
    subplot(2+do_direct,1,1);
    plot(Z{i},{'x','y'},'b',0.2)
    title('Memory Approach')
    subplot(2+do_direct,1,2);
    plot(Z{i},{'x','y'},'k',0.05)
    plot(liftZ_newdata,{sprintf('x_%i',i),sprintf('y_%i',i)},'r',0.2)
    title('Memory with new data at k=3')
    if do_direct
        subplot(3,1,3);
        plot([1,0,0;0,1,0]*Z_direct{i},'b',0.2);
        title('non-Memory Approach')
    end
end
subplot(2+do_direct,1,1);
ax1 = gca();
subplot(2+do_direct,1,2);
ax2 = gca();
ax2.XLim = ax1.XLim;
ax2.YLim = ax1.YLim;
if do_direct
    subplot(2+do_direct,1,3);
    ax3 = gca();
    ax3.XLim = ax1.XLim;
    ax3.YLim = ax1.YLim;
end

%% Helper Functions
function [Z] = evolveUnicycleMemZono(sinZono,cosZono,v,U,Z0,N)
    Z = {memZono(Z0,sprintf('X_1'))};
    % Give the state the correct dimension labels
    Z{1}.dimKeys = {'x','y','theta'};
    
    for k = 2:N
        %===== These need to be inside the loop so that each usage of sin/cos/U is memory independent
        sinMZ = memZono(sinZono,sprintf('sin_k%i',k));
        % sin theta updates the y state
        sinMZ.dimKeys = {'theta','y'};
    
        cosMZ = memZono(cosZono,sprintf('cos_k%i',k));
        % cos theta updates the x state
        cosMZ.dimKeys = {'theta','x'};
    
        UMZ = memZono(U,sprintf('U_%i',k-1));
        UMZ.dimKeys = {'theta'};
        %=====
    
        % does X_{k+1}  = X_k
        %      Y_{k+1}  = Y_k
        %      TH_{k+1} = TH_K + U
        Z{k} = combine(Z{k-1},UMZ);
    
        % does X_{k+1} = ... + v*cos(theta)
        % only use Z{k-1}({'theta'}) otherwise the x and y dims would get intersected
        dX = v*merge(cosMZ,Z{k-1}({'theta'}),sprintf('cos_x_%i',k));
        Z{k} = combine(Z{k}, dX({'x'})); % only add the x (the function value, not the theta input domain of cosine)
    
        % does Y_{k+1} = .. + v*sin(theta)
        dY = v*merge(sinMZ,Z{k-1}({'theta'}),sprintf('sin_y_%i',k));
        Z{k} = combine(Z{k}, dY({'y'}));
    end
end

function [X] = evolveUnicycle(sinZono,cosZono,v,U,X0,N)
    X = {X0};
    
    for k = 2:N
        XYTH_prev = X{k-1};
        % set of xs at previous time step
        X_prev = [1,0,0]*XYTH_prev;
        % set of ys at previous time step
        Y_prev = [0,1,0]*XYTH_prev;
        % set of thetas at previous time step
        TH_prev = [0,0,1]*XYTH_prev;
    
        % use [0,1] to only keep the sin/cos function output
        X_next = X_prev + v*[0,1]*and(cosZono, TH_prev, [1,0]);
        Y_next = Y_prev + v*[0,1]*and(sinZono, TH_prev, [1,0]);
        TH_next = TH_prev + U;
    
        X_k = cartProd( cartProd(X_next,Y_next),TH_next);
        X{k} = X_k;
    
        %plot([1,0,0;0,1,0]*X_k,'b',0.3)
    end
    hold off;
end

function cosx = makeSinX(n,bd)

    th = linspace(bd(1), bd(2), 2*n+1);
    xi = th;
    yi = sin(th);
    
    V = [xi; yi];
    nc = 2*n+1;
    nb = nc-1;
    naux = nc;
    c = zeros(2,1);
    Gc = V;
    
    M = zeros(nc,nb);
    for i = 1:nb
        M(i,i) = 1;
        M(i+1,i) = 1;
    end
    
    c = c+0.5*Gc*ones(nc,1);
    Gc = 0.5*Gc;
    Gc = [Gc,zeros(2,naux)];
    Gb = zeros(2,nb);
    
    Ac = [0.5*ones(1,nc) zeros(1,naux)]; %(1)
    Ab = zeros(1,nb);
    b = [1-0.5*nc];
    
    Ac = [Ac; 0.5*eye(nc) 0.5*eye(nc)]; %(2)
    Ab = [Ab; -0.5*M];
    b = [b; 0.5*M*ones(nb,1)-ones(nc,1)];
    
    Ac = [Ac; zeros(1,nc+naux)]; %(3)
    Ab = [Ab; 0.5*ones(1,nb)];
    b = [b; 1-0.5*nb];
    
    cosx = hybZono(Gc,Gb,c,Ac,Ab,b);
    
end

function cosx = makeCosX(n,bd)

    th = linspace(bd(1), bd(2), 2*n+1);
    xi = th;
    yi = cos(th);
    
    V = [xi; yi];
    nc = 2*n+1;
    nb = nc-1;
    naux = nc;
    c = zeros(2,1);
    Gc = V;
    
    M = zeros(nc,nb);
    for i = 1:nb
        M(i,i) = 1;
        M(i+1,i) = 1;
    end
    
    c = c+0.5*Gc*ones(nc,1);
    Gc = 0.5*Gc;
    Gc = [Gc,zeros(2,naux)];
    Gb = zeros(2,nb);
    
    Ac = [0.5*ones(1,nc) zeros(1,naux)]; %(1)
    Ab = zeros(1,nb);
    b = [1-0.5*nc];
    
    Ac = [Ac; 0.5*eye(nc) 0.5*eye(nc)]; %(2)
    Ab = [Ab; -0.5*M];
    b = [b; 0.5*M*ones(nb,1)-ones(nc,1)];
    
    Ac = [Ac; zeros(1,nc+naux)]; %(3)
    Ab = [Ab; 0.5*ones(1,nb)];
    b = [b; 1-0.5*nb];
    
    cosx = hybZono(Gc,Gb,c,Ac,Ab,b);
    
end