% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example:
%   Forwared reachability with a single guard condition and two equilibrium
%   points.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Parameters

N = 15;     % Number of discrete time steps
bounds = 3; % Bounds on state space for analysis (assuming [-bounds bounds] in each dimension)

% System dynamics
Ad1 = [0.75,-0.25;0.25,0.75];
Bd1 = [0.25;-0.25];
Ad2 = [0.75,0.25;-0.25,0.75];
Bd2 = [-0.25;-0.25];

% Initialize plotting
figure; hold on;
axis([ -2 2 -1 3 ]);
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
set(gca,'fontsize',18,'fontname','times new roman')
plot(-1,0, 'k.','markersize',25) % Left equilibrium
plot(1,0, 'k.','markersize',25)  % Right equilibrium
p = plot([ 0 0 ] , [ -bounds bounds ], 'g--','linewidth',2); % Guard
colors = interp1([1;N+1],[0 0 1;1 0 0],1:1:N+1);

% Find a good initial condition set
cx = [ -0.201 ; 1.395 ];	% Crosses the guard once
Gx = 0.2*eye(2);
Z0 = zono(Gx,cx);	
% Propagate backwards two steps
Ainv = sysD2.A^(-1);
Z0 = Ainv*(Z0 + -1*sysD2.B);
Z0 = Ainv*(Z0 + -1*sysD2.B);
plot(Z0,colors(1,:),1) % Plot initial condition set

% No input set
U = [];

% Generated by Hysdel
A = zeros(2);
B_u = [];
B_w = [1 0 0; 0 1 0];
B_aff = [0;0];
E_x = [1 0; -1 0; 0.75 -0.25; -0.75 0.25; 0.75 0.25; -0.75 -0.25; 0.25 0.75; -0.25 -0.75; -0.25 0.75; 0.25 -0.75];
E_u = [];
E_w = [0 0 -5; 0 0 5; -1 0 10.5; 1 0 9.5; -1 0 -9.5; 1 0 -10.5; 0 -1 10; 0 1 10; 0 -1 -10; 0 1 -10];
E_aff = [0; 5; 10.25; 9.75; 0.25; -0.25; 10.25; 9.75; 0.25; -0.25];

wUB = [5.25;4.74;1];
wLB = [-5.25;-5.25;0];

W = hybZono([diag((wUB(1:2)-wLB(1:2))/2);zeros(1,2)],[zeros(2,1);diag((wUB(3)-wLB(3))/2)],(wUB+wLB)/2,[],[],[]);

% Compute and plot reachable sets
Z = Z0;
for i = 1:N
    Z = stepMLD(Z,U,W,A,B_u,B_w,B_aff,E_x,E_u,E_w,E_aff);
    plot(Z,colors(i+1,:),1) % Plot reachable set
end

% Additional plot information
legend(p,{'Guard'},'interpreter','latex')

%% Testing
% figure; hold on
% X0 = zono(bounds*eye(2),zeros(2,1));
% plot(X0,'k',0.1)
% X0a = halfspaceIntersection(X0,[1 0],0);
% X0b = halfspaceIntersection(X0,[-1 0],0);
% plot(X0a,'m',0.1)
% plot(X0b,'g',0.1)
% X1a = Ad1*X0a + Bd1;
% plot(X1a,'m',0.1)
% X1b = Ad2*X0b + Bd2;
% plot(X1b,'g',0.1)
% plot(union(X1a,X1b),'k',0.5)