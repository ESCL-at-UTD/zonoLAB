clear; close all;

x_max = 1.5*pi;
x_min = -x_max;

%% Sine Function as Hybrid Zonotope
n_points = 11; % Number of points (should be an odd number)
x = linspace(x_min,x_max,n_points);
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
sinSOS = makeSinX(n_points,[x_min x_max]);
sinZono = sinSOS+sinErr;

axislim = [x_min,x_max,-1.1,1.1];
figure(1);
set(gcf,'Position',[100 100 500 1200])
subplot(4,1,1)
fplot(@(x) sin(x),[x_min x_max],'-','color','g','LineWidth',2)
plot(sinZono,'b',1)
hold on
plot(sinSOS,plotOptions('EdgeColor','r','FaceColor','r','LineWidth',2))
fplot(@(x) sin(x),[x_min x_max],'-','color','g','LineWidth',2)
grid on;
axis(axislim)
legend('$\sin(x)$',...
        '$\bar{\mathcal{Z}}_{\sin(x)}$',...
        '$\mathcal{Z}_{SOS2}$',...
        'Interpreter','latex','location','best')
hold off;
title('Approximation of Sine','Interpreter','latex')
%set(gca,'FontSize',24,'FontName','times')

%% Using the Hybrid Zonotope approximation of Sine as a function

subplot(4,1,2)
% Define a domain over x of [-pi/2,pi/2]
X1 = zono([pi/2],0);
% select the part of the hybZono sine function set that intersects the
% domain X1. The generalized intersection uses [1,0] to intersect only in
% the first dimension of the sine function set.
sin_X1 = and(sinZono,X1,[1,0]);
plot(X1);
plot(sin_X1,'m',0.7)
title('$\cap$ with $x=[-\pi/2,\pi/2]$','Interpreter','latex')
axis(axislim)

subplot(4,1,3)
% Define a range over y of [0,0.5]
Y1 = zono([0.25],0.25);
% Define this range so we can plot it for visualization.
% note that this definition explicitly sets x=0 (which is not how we are using it for the intersection below)
Y1_forplot = zono([0;0.25],[0;0.25]); 
% select the part of the hybZono sine function set that intersects the
% range Y1. The generalized intersection uses [0,1] to intersect only in
% the second dimension of the sine function set.
sin_Y1 = and(sinZono,Y1,[0,1]);
plot(Y1_forplot);
plot(sin_Y1,'y',0.7)
title('$\cap$ with $y=[0,0.5]$','Interpreter','latex')
axis(axislim)

subplot(4,1,4)
% Define a disconnected range over y of [-0.75,-0.25] and [0.25,0.75]
Y2 = hybZono([0.25],[0.5],0,[],[],[]);
Y2_forplot = hybZono([0;0.25],[0;0.5],[0;0],[],[],[]);
sin_Y2 = and(sinZono,Y2,[0,1]);
plot(Y2_forplot)
plot(sin_Y2,'g',0.7)
title('$\cap$ with $y\in\{[-0.75,-0.25]\cup[0.25,0.75]$','Interpreter','latex')
axis(axislim)


function sinx = makeSinX(n,bd)

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
    
    sinx = hybZono(Gc,Gb,c,Ac,Ab,b);
    
end