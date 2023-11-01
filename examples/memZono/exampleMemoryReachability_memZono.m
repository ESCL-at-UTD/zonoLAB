clear; clc; clf;

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

% Initial and Final Sets (really need to get constructor better)
X_0 = memZono(zono(diag([1,2]),zeros(2,1)),'x_0');
X_F = memZono(zono(0.25*eye(2),ones(2,1)),'x_f');

% Input Set
U_nom = zono(1.5,0);
U_0 = memZono(U_nom,'u_0');
for k = 1:N
    U_{k} = memZono(U_nom,sprintf('u_%d',k));
end

%% Reachability Calculation
% Initialization
Xall = [X_0; U_0]; %<--- cart prod but extendable to many of them
X_{1} = A*X_0 + B*U_0;
for k = 1:N
    % Name and save/label
    X_{k}.dimKeys = sprintf('x_%d',k);
    Xall = [Xall; X_{k}; U_{k}];
    % Next Time-step
    X_{k+1} = A*X_{k} + B*U_{k};
end
% Intersection
X_F.dimKeys = sprintf('x_%d',k);
Xall = Xall & X_F 

%% 



% %% Simulation Settings
% 
% 
% % Time
% T = 2*pi;  % this will get slower the longer you run it because the X zonotope will get more and more complex!
% dt = 0.25;
% N = round(T/dt)+1;
% time = linspace(0,T,N);
% 
% %%  System definition
% % Spring-Mass-Damper:
% A = [0,1;-2,-1];
% B = [0;1];
% 
% F = expm(A*dt); % Discrete-time dynamics
% [n,m] = size(B);
% 
% tau = linspace(0,dt,1001); % integration steps for calculating G
% expmAtauB = zeros(n,m,length(tau));
% for i=1:length(tau)
%     expmAtauB(:,:,i) = expm(A*tau(i))*B;
% end
% G = trapz(tau, expmAtauB, 3); % numerically integrate (along dimension 3)
% 
% eta = 1; % noise level - noise will be pulled from [-eta,eta]
% 
% % select estimator gain to be stable - for discrete systems it means
% % eigenvalues between (-1,1):
% C = eye(2); % we assume we'll measure both states
% L = place(F',C',0.8*[1,1]); %%%%% You can try different eigenvalues here between (-1,1)
% 
% % select controller gain to be stable
% K = -place(F,G,[0.5;0.4]);
% 
% %% Test
% Ns = 3;
% Unom = conZono(1.5,0,[],[]); % just a big set for input
% X={}; X{1} = conZono(diag([1,2]),zeros(2,1),[],[]);
% XF = conZono(0.25*diag([1,1]),ones(2,1),[],[]);
% 
% XX = X{1};
% rx = {}; rx{1} = 1:n; % use for indexes that correspond to states
% ru = {}; % use for indexes that correspond to inputs
% ri = n; % last used index
% for k = 2:Ns
%     % Compute input set from Hybrid Zonotope verison Neural Net
%     U{k-1} = Unom;
%     % Compute AX + BU
%     X{k} = plus( F*X{k-1} , G*U{k-1} );
%     %XU = X{k};
%     XU = extend_zonotope(X{k},U{k-1});
% %     % enforcing the intersection along the way
% %     if k == Ns
% %         R = zeros(n,length(XU.c));
% %         R(:,1:n) = eye(n);
% %         XU = generalizedIntersection(XU,XF,R);
% %     end
%     rx{k} = ri+1:ri+n; ri = ri+n;
%     ru{k-1} = ri+1:ri+m; ri = ri+m;
% 
%     % Now stack X_k and X_{k+1} so that we keep factors across time
%     XX = extend_zonotope(XX,XU);
% end 
% % if we want to do the intersection at the end
% R = zeros(n,length(XX.c));
% R(:,rx{3}) = eye(n);
% XX = and(XX,XF,R);
% 
% %% Plot Sets
% figure(1)
% clf
% 
% subplot(1,2,1)
% plot(XF, 'g', 1);
% drawnow;
% hold on;
% for j = 1:Ns  % = 1,3,5,7,9,11 to select the different time steps
%     clr='k'; 
%     if mod(j,2) == 0
%         clr = 'b';
%     elseif j == 1
%         clr = 'k'; 
%     else
%         clr = 'r';
%     end
%     R = zeros(n,length(XX.c));
%     R(:,rx{j}) = eye(n);
%     plot( R*XX, clr, 0.6)
%     hold on;
%     plot( X{j}, clr, 0.2)
%     drawnow;
% 
% end
% hold off;
% 
% axis equal;
% xlim([-3 3]);
% ylim([-3 3]);
% xlabel('$x_1$','Interpreter','latex');
% ylabel('$x_2$','Interpreter','latex');
% 
% subplot(1,2,2)
% R = zeros(2,length(XX.c));
% R(1,ru{1}) = 1;
% R(2,ru{2}) = 1;
% plot( R*XX, 'b', 0.6)
% hold on;
% plot( stack_zonotopes(U{1},U{2}) , 'b', 0.2)
% hold off;
% drawnow;
% % for j = 1:Ns-1
% %     R = zeros(m,length(XX.c));
% %     R(:,ru{j}) = eye(m);
% %     Ui = R*XX;
% %     Ui.c = [j;Ui.c];
% %     Ui.Gc = [j*zeros(1,size(Ui.Gc,2));Ui.Gc];
% %     Ui.Gb = [j*zeros(1,size(Ui.Gb,2));Ui.Gb];
% %     plot( Ui, 'b', 0.6)
% %     drawnow;
% %     hold on;
% % end
% hold off;
% axis equal;
% xlim([-2 2]);
% ylim([-2 2]);
% 
% xlabel('$u(1)$','Interpreter','latex');
% ylabel('$u(2)$','Interpreter','latex');
% 
% %% Functions
% 
% function [Z] = extend_zonotope(X,Y)
% % Extending X with Y. 
% nx = size(X.G,1);
% ny = size(Y.G,1);
% ncx = size(X.A,1);
% ncy = size(Y.A,1);
% r_ngc_diff = 0;
% l_ngc_diff = 0;
% if size(Y.G,2) >= size(X.G,2)
%     %If X has f factors, then the first f factors of Y are the same as X 
%     r_ngc_diff = size(Y.G,2) - size(X.G,2); % num new cont. factors
% elseif size(X.G,2) > size(Y.G,2)
%     %If Y has f factors, then the last f factors of X are the same as Y 
%     l_ngc_diff = size(X.G,2) - size(Y.G,2); % num new cont. factors
% else
%     error(sprintf('Ambiguous extension. One set should have more generators than the other (both binary and continuous).')) 
% end
% 
% Z = conZono( [X.G , zeros(nx,r_ngc_diff) ; zeros(ny,l_ngc_diff) , Y.G, ], ...
%              [X.c ; Y.c], ...
%              [X.A , zeros(ncx,r_ngc_diff) ; zeros(ncy,l_ngc_diff) , Y.A], ...
%              [X.b ; Y.b]);
% end
% 
% function [Z] = stack_zonotopes(X,Y)
% % Stacking two independent zonotopes. 
% nx = size(X.G,1);
% ny = size(Y.G,1);
% ncx = size(X.A,1);
% ncy = size(Y.A,1);
% 
% Z = conZono( [X.G , zeros(nx,size(Y.G,2)) ; zeros(ny,size(X.G,2)) , Y.G, ], ...
%              [X.c ; Y.c], ...
%              [X.A , zeros(ncx,size(Y.G,2)) ; zeros(ncy,size(X.G,2)) , Y.A], ...
%              [X.b ; Y.b]);
% end