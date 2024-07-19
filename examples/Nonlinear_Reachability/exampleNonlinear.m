% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example:
%   Closed-loop reachability of a nonlinear system with a saturated LQR
%   controller.
%   System definition:
%   x1[k+1] = x1[k] + dt*x2[k]
%   x2[k+1] = x2[k] + dt*(10*sin(x1[k]) + u[k])
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 
tAll = tic;
%% Parameters
dt = 0.1;   % Time step size
n = 2;      % Number of states
m = 1;      % Number of inputs
u_max = 20; % Maximum input
u_min = -u_max;
x_max = [4;8]; % Maximum state
x_min = -x_max;
% Linearized model about x = 0
A = [1 dt; 10*dt 1];
B = [0; dt];
% LQR design
Q = eye(n);
R = eye(m);
[K,~,~] = dlqr(A,B,Q,R);
K = -K; % Assuming u = K*x
% State domain 
X = zono(diag(x_max),zeros(n,1));

%% Saturated LQR
% LQR set-based map
UX_LQR = [eye(2);K]*X;
% Saturation 
M = sum(abs(K*X.G)); % Maximum value of u = K*x for x_min <= x <= x_max
cL = (M+u_max)/2;
gL = (M-u_max)/2;
Prop = zono(u_max*[1 1]',zeros(2,1));
Sat = hybZono([gL;0],[cL;u_max],zeros(2,1),[],[],[]);
PropSat = union(Prop,Sat);
% Saturated LQR
UX_LQR_SAT = cartProd(UX_LQR,PropSat);
UX_LQR_SAT = coupleStates(UX_LQR_SAT,[3,4]);
UX_LQR_SAT = [1 0 0 0 0; 0 1 0 0 0; 0 0 0 0 1]*UX_LQR_SAT;
% Plot control law
figure;
plot(UX_LQR_SAT,'r',1)
xlabel('$x_{1,k}$','Interpreter','latex')
ylabel('$x_{2,k}$','Interpreter','latex')
zlabel('$u(x_k)$','Interpreter','latex')
grid on
set(gca,'FontSize',24,'FontName','times')
axis([-4 4 -8 8 -20 20]);
ax = gca;
ax.XTick = [-4 0 4];
ax.YTick = [-8 0 8];
set(gcf,'Position',[100 100 400 350]);
view(70,30)

exportgraphics(gcf,'nonlin_LQR.pdf','ContentType','vector')
%% Sine Function as Hybrid Zonotope
n_points = 21; % Number of points (should be an odd number)
x = linspace(x_min(1),x_max(1),n_points);
y = sin(x);

errSOS = zeros(n_points-1,2);
for i = 1:n_points-1
    x0 = x(i);
    x1 = x(i+1);
    a = (sin(x1)-sin(x0))/(x1-x0);
    b = (x1*sin(x0) - x0*sin(x1))/(x1-x0);
    if 2*pi-acos(a) >= x0 && 2*pi-acos(a) <= x1
        errSOS(i,1) = a*(2*pi-acos(a)) + b + sin(acos(a));
    end
    if acos(a) >= x0 && acos(a) <= x1
       errSOS(i,2) = -a*acos(a) - b + sin(acos(a));
    end
end
errSOS_max = max(abs(errSOS(:)));

sinErr = zono([0;errSOS_max],zeros(2,1));
sinSOS = makeSinX(n_points,[x_min(1) x_max(1)]);
sinSUS = sinSOS+sinErr;


figure;
fplot(@(x) sin(x),[x_min(1) x_max(1)],'-','color','g','LineWidth',2)
plot(sinSUS,'b',1)
hold on
plot(sinSOS,plotOptions('EdgeColor','r','FaceColor','r','LineWidth',2))
fplot(@(x) sin(x),[x_min(1) x_max(1)],'-','color','g','LineWidth',2)
grid on;
TMP = [1.1 2.1 .9 1.05];
plot([TMP(1) TMP(2)],TMP(3)*ones(1,2),'k-','LineWidth',3)
plot([TMP(1) TMP(2)],TMP(4)*ones(1,2),'k-','LineWidth',3)
plot(TMP(1)*ones(1,2),[TMP(3) TMP(4)],'k-','LineWidth',3)
plot(TMP(2)*ones(1,2),[TMP(3) TMP(4)],'k-','LineWidth',3)

xlabel('$x$','Interpreter','latex')
ylabel('$y$','Interpreter','latex')
set(gca,'FontSize',24,'FontName','times')
set(gcf,'Position',[100 100 900 400])
legend('$\sin(x)$','$\bar{\mathcal{Z}}_{\sin(x)}$','$\mathcal{Z}_{SOS2}$','Interpreter','latex','location','best')
hold off

axes('Position',[0.55 0.35 0.35 0.35])
fplot(@(x) sin(x),[x_min(1) x_max(1)],'-','color','g','LineWidth',2)
plot(sinSUS,'b',1)
hold on
plot(sinSOS,plotOptions('EdgeColor','r','FaceColor','r','LineWidth',2))
fplot(@(x) sin(x),[x_min(1) x_max(1)],'-','color','g','LineWidth',2)
grid on;
axis([1.1 2.1 .9 1.05])
plot([TMP(1) TMP(2)],TMP(3)*ones(1,2),'k-','LineWidth',3)
plot([TMP(1) TMP(2)],TMP(4)*ones(1,2),'k-','LineWidth',3)
plot(TMP(1)*ones(1,2),[TMP(3) TMP(4)],'k-','LineWidth',3)
plot(TMP(2)*ones(1,2),[TMP(3) TMP(4)],'k-','LineWidth',3)
set(gca,'FontSize',24,'FontName','times')

% exportgraphics(gcf,'nonlin_sin_function.pdf','ContentType','vector')
%% Open-loop State-Update Set 
Phi_OL = zono([4 0 0 0; 0 3*x_max(2) 0 0; 0 0 u_max 0; 0 0 0 0; 0 0 0 10],zeros(5,1));
Phi_OL = cartProd(Phi_OL,diag([1 10])*sinSUS);
Phi_OL = coupleStates(Phi_OL,[1 6]);
Phi_OL = coupleStates(Phi_OL,[5 7]);
Phi_OL = [eye(5) zeros(5,2)]*Phi_OL; %[x1_k x2_k T_k x1_k+1 sin(x1_k)]

Phi_OL = [  1 0 0 0 0;
            0 1 0 0 0;
            0 0 1 0 0;
            1 dt 0 0 0;
            0 1 dt 0 dt]*Phi_OL;

%% Closed-loop State-Update Set
Phi_CL = [eye(2) zeros(2,3); zeros(2,3) eye(2)]*and(Phi_OL,UX_LQR_SAT,[eye(3) zeros(3,2)]);

figure
PhiNL_CL_124 = projection(Phi_CL,[1 2 4]);
plot(PhiNL_CL_124,'b',1)
xlabel('$x_{1,k}$','Interpreter','latex');
ylabel('$x_{2,k}$','Interpreter','latex');
zlabel('$x_{2,k+1}$','Interpreter','latex');
grid on;
set(gca,'FontSize',16,'FontName','times')
axis([-4 4 -8 8 -11 11]);
ax = gca;
ax.XTick = [-4 0 4];
ax.YTick = [-8 0 8];
set(gcf,'Position',[100 100 500 400]);
view(70,30)

%% Closed-Loop Reachability Analysis

n_Sim = 12; % Number of discrete time steps to compute

Xreach{1} = zono(diag([pi 0.1]),[0;0]); % Initial condition set
complexity = [Xreach{1}.nG 0 0 0 0]; % [nGc nGb nC n_sets tCalc]
ZeroEye = [zeros(2) eye(2)];
EyeZero = [eye(2) zeros(2)];
for i = 1:n_Sim
    Xreach{i+1} = ZeroEye*and(Phi_CL,Xreach{i},EyeZero);
%     leaves = getLeaves(Xreach{i+1},solverOptions);
%     nLeaves = size(leaves,2);
%     complexity(i+1,:) = [Xreach{i+1}.nGc Xreach{i+1}.nGb Xreach{i+1}.nC nLeaves 0];
    i
end

%% Plot Reachable Sets
tPlot = tic;
n_plot =n_Sim;
% Colors for sets
blu = [0 0 1]; % Blue
red = [1 0 0]; % Red
Cx = 1:1:n_plot+1; % Values to interpolate on 
CxOk = [1;n_plot+1];
colors = interp1(CxOk,[blu;red],Cx);
% Color bar
tickLab{1} = ['$\mathcal{X}_{',num2str(0),'}$'];


figure
for i = 1:n_plot+1
    tStart = tic;
    plot(Xreach{i},plotOptions('EdgeColor',colors(i,:),'FaceColor',colors(i,:),'LineWidth',0.01))
    drawnow;
    complexity(i,5) = toc(tStart);
    if i == 1
        tickLab{i} = ['$\mathcal{X}_{',num2str(0),'}$'];
    else
        tickLab{i} = ['$\mathcal{R}_{',num2str(i-1),'}$'];
    end
end
colormap(colors);
cbar = colorbar;
cbar.TicksMode = 'manual';
val = 1/(((n_plot+1)));
cbar.Ticks = val/2+[0:val:1-val];
cbar.TickLabelInterpreter = 'latex';
cbar.TickLabels = tickLab;
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
xlim([-4 4])
ylim([-7 7])
grid on
fontsize = 24;
set(gca,'FontSize',fontsize,'fontname','times new roman')
toc(tPlot)
toc(tAll)
% exportgraphics(gcf,'nonlin_CL_reach.pdf','ContentType','vector')

%% Local Functions
function [out] = coupleStates(in,e)

Ac_add = in.Gc(e(:,2),:)-in.Gc(e(:,1),:);
Ab_add = in.Gb(e(:,2),:)-in.Gb(e(:,1),:);
b_add = in.c(e(:,1))-in.c(e(:,2));

out = in;
out.Ac = [out.Ac; Ac_add];
out.Ab = [out.Ab; Ab_add];
out.b = [out.b; b_add];

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