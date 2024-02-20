% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example:
%   Set Containment Reachability Optimization with ConZono
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

clear all; close all; clc;
%% LTI System

A = [0.84 0.23; -0.47 0.12];
B = [0.07 0.3; 0.23 0.1];
Nx = size(A,1);
Nu = size(B,2);

% Define Dangerous State Region
uBound = [15 10];
lBound = [-15 0];
d = [-0.1 -1];
b = -5;
D_HRep = Polyhedron('H', [d, b], 'lb', lBound, 'ub', uBound);

% Define Safe State Region
S.G = [0 -200; 30 20];
S.c = [0 -25]';
S.A = [0 0; 0 0; 0 0; 0 0];
S.b = [0 0 0 0]';
S = conZono(S.G, S.c, S.A, S.b);

% Define Natural Input set
U_n.G = diag([8,6]); 
U_n.c = [-1 -1]';
U_n = zono(U_n.G, U_n.c);

% Define Designed Input set
U_d.G = [1 0 1; 0 1 -1];
phi = diag([3.76, 2.06, 3.94]);
U_d.c = [-1.3 -1]';
U_d = zono(U_d.G*phi, U_d.c);

% Define Operational States
O1 = [-8.044 -7.6812 -7.1369 -6.5564 -5.8308 -5.0507 -5.1233 -5.7764 -6.4839 -7.4453];
O2 = [3.4254 4.3869 4.4594 4.2962 4.0785 3.8608 3.6431 3.3710 3.3165 3.2440];



%% Operational Set O
O_HRep = Polyhedron('V', [O1' O2']);
H_O = O_HRep.H(:,1:end-1);
k_O = O_HRep.H(:,end);

O = hPoly2conZono([H_O, k_O]);

%% Reachability Set R_N
N = 10;

%Initial Conditions
I = eye(size(A,1));
X0.G = (I-A)\B*U_d.G*phi;
X0.c = (I-A)\B*U_d.c;
X0 = zono(X0.G, X0.c);

R_N = X0;

%Compute Reachable Set using Zonotopes
for i = 1:N
    R_N = A*R_N + B*U_n;
end

R_N = conZono(R_N);

%% Designed set R_D

%Initial Conditions
X0.G = zeros(2,0);
X0.c = zeros(2,1);
X0 = zono(X0.G, X0.c);

R_D = X0;

%Compute Reachable Set using Zonotopes
for i = 1:N
    R_D = A*R_D + B*U_d;
end

R_D = conZono(R_D);

%% Optimization: Safe Set / Reachability / Operational set containments

%Create Constraints for Operational set O within Desired Reachable set R_D
reachable = optimproblem;
reachable.Objective = 0;
reachable = enforceSetContain(O, R_D, reachable);

%Create Constraints for Desired Reachable Set R_D within Safe set S
reachable = enforceSetContain(R_D, S, reachable,'2');

%Solve Optimization problem
options = optimoptions('linprog','Display','off');
[sol, fval, exitflag] = solve(reachable, 'Options',options);

if exitflag == 1
    disp('    RESULT: Set is contained.')
else
    disp('    RESULT: Set is not contained.')
end

%Plot
figure();
plot(D_HRep,'color','r','alpha', 0.1, 'LineStyle', 'none')
plot(S, 'g', 0.1);
plot(R_D,'k',0.4);
hold on;
plot(O,'b',0.5);
xlim([-15 15]);
ylim([-10 10]);
xlabel('$x_1$','interpreter', 'Latex');
ylabel('$x_2$','interpreter', 'Latex');
legend('D', 'S', 'R_D', 'O')
grid on;
box on;