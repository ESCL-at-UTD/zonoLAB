% mpt_demo_lti4
A = [1 1; 0 1];
B = [1; 0.5];
C = [1 0];
D = 0;
model = LTISystem('A', A, 'B', B, 'C', C, 'D', D);

% Next we define an MPC controller:
horizon = 5;
ctrl = MPCController(model, horizon);

% Specify the MPC problem (constraints and penalties):
ctrl.model.x.min = [-5; -5];
ctrl.model.x.max = [5; 5];
ctrl.model.u.min = -1;
ctrl.model.u.max = 1;
% Use quadratic cost function
ctrl.model.x.penalty = QuadFunction(eye(model.nx));
ctrl.model.u.penalty = QuadFunction(eye(model.nu));

% Finally, we can convert the controller to an explicit form:
disp('Generating the explicit solution:');
tStart = tic;
expctrl = ctrl.toExplicit();
toc(tStart)

% Compare optimal solutions
% x0 = [-4; 0];
% disp('Optimal control input obtained by evaluating the explicit solution:');
% Uexp = expctrl.evaluate(x0)
% 
% disp('Optimal control input obtained by solving the optimization problem on-line:');
% Uonl = ctrl.evaluate(x0)

% plot the explicit optimizer
% figure
% expctrl.feedback.fplot();
% title('PWA representation of the optimal control input')
% 
% % plot the value function
% figure
% expctrl.cost.fplot();
% title('Explicit cost function');

% plot the regions
figure
expctrl.partition.plot()
title('Regions of the polyhedral partition');

%% zonoLab Implementation
%   min_U 1/2 U'*H*U + x'*F*U + Cf*U + x'*Y*x + Cx*x + Cc
%    s.t. G*U <= W + E*x
Matrices = getMatrices(ctrl);

% dimensions of the mp-QP
nz = size(Matrices.H,2);
nx = size(Matrices.E,2);
nh = size(Matrices.G,1);

% represent the mp-QP as a hybrid zonotope

tStart = tic;

% define parameter set over all of the bounded parameters
Gx = diag((ctrl.model.x.max-ctrl.model.x.min)/2);
cx = (ctrl.model.x.max+ctrl.model.x.min)/2;
X = zono(Gx,cx);
% figure;
% plot(X)
X = hybZono(X);

% find objective function and constraints in the form
% min_U 1/2U'*Q*U + x'*P*U + q'*U
%	s.t. HU + Sx <= f

Q = Matrices.H;
P = Matrices.F;
q = Matrices.Cf';
H = full(Matrices.G);
S = full(-Matrices.E);
f = full(Matrices.W);

% Find min H-rep
feasSet = Polyhedron('H',[H S f]);
feasSet.minHRep;
nh
H = feasSet.H(:,1:nz);
S = feasSet.H(:,1+nz:end-1);
f = feasSet.H(:,end);
nh = size(H,1)

% now convert to standard form
% min_z 1/2z'*Q*z + q'*z
%	s.t. Hz + Sx <= f
%	where z = U - inv(Q)*P*x 
S = S - H*inv(Q)*P';

% check if the problem is strictly convex
[~,p1] = chol(Q);
if p1
	error('The objective function is not positive definite..')
end

% Find over approximative zonotope Zb > Z*
% z = U + Px
% s.t. Hz + Sx <= f
% x in X -> over approximate as x in parambounds

f_max = f - S*X.c + sum(abs(S*X.Gc),2);

Z = conZono([H f_max]);
% figure;plot(projection(Z,[2 3 5]),'b',0.1)
Zb = boundingBox(Z);
Gz = Zb.G;
cz = Zb.c;

%
% formulate the hybrid zonotopes Z* and Xfeas

% M = 1e3*ones(nh,1);
muMaxT = 1e5*ones(nh,1);

% Compute bounds
ZX = conZono([H S f]);
ZX = and(ZX,X,[zeros(nx,nz) eye(nx)]);

ZX_mapped = [H S]*ZX + zono(-f);
ZX_mappedZono = zono(ZX_mapped.Gc,ZX_mapped.c); % Conver to zonotope overapprox set and reduce computation time for boundingBox
ZX_mapped_box = boundingBox(ZX_mapped);
M = abs(ZX_mapped_box.c - ZX_mapped_box.G*ones(nh,1));


% for now take bounds on big-M of mu as max value.. don't want to accidentally trim feasible space just yet
% muMaxTt = muMaxT;
% muMaxT = 10*max(abs(muMaxT))*ones(nh,1);	% factor of 10

% transform objective and constraints to factor space
Qb = Gz'*Q*Gz;
qb = (cz'*Q*Gz + q'*Gz)';
Hb = H*Gz;
fb = f - H*cz;
Hmu = Hb'*diag(muMaxT/2);
Ss = [ S ; -1*S ; zeros(nh,nx) ];
Hs = [ Hb ; -1*Hb ; zeros(nh,nz) ];
Ts = [ zeros(2*nh,nh) ; eye(nh) ];
Cs = [ zeros(nh) ; diag(M/2) ; -1*eye(nh) ];
fs = [ fb ; -fb + M/2 ; zeros(nh,1) ];

% now build the hybrid zonotope
% left limits on the inequality constraints
alpha = -sum(abs([ Hb S*[X.Gc X.Gb] ]),2);
beta = -sum(abs([ Hb S*[X.Gc X.Gb] ]),2) - M/2;
fbLeft = [ alpha ; beta ; -2*ones(nh,1) ];
fbRight = fs - Ss*X.c;
Gs = (fbLeft - fbRight)./2;
Gs = diag(Gs);
cs = (fbRight + fbLeft)./2;

% constraints
Acs1 = [ X.Ac zeros(X.nC,nz+4*nh) ];
Abs1 = [ X.Ab zeros(X.nC,nh) ];
bs1 = X.b;
Acs2 = [ zeros(nz,X.nGc) Qb Hmu zeros(nz,3*nh) ];
Abs2 = zeros(nz,X.nGb+nh);
bs2 = -qb - Hmu*ones(nh,1);
Acs3 = [ Ss*X.Gc Hs Ts Gs ];
if ~isempty(X.Gb)
	Abs3 = [ Ss*X.Gb Cs ];
else
	Abs3 = Cs;
end
bs3 = cs;

% all together
Acs = [ Acs1 ; Acs2 ; Acs3 ];
Abs = [ Abs1 ; Abs2 ; Abs3 ];
bs = [ bs1 ; bs2 ; bs3 ];

Xfeas = X;
Xfeas.Gc = [ X.Gc zeros(nx,nz+4*nh) ];
Xfeas.Gb = [ X.Gb zeros(nx,nh) ];
Xfeas.Ac = Acs;
Xfeas.Ab = Abs;
Xfeas.b = bs;

toc(tStart)

figure;
plot(Xfeas,'b',0.1)
%%
R = rref([Xfeas.Ac Xfeas.Ab Xfeas.b]);
Xfeas.Ac = R(:,1:Xfeas.nGc);
Xfeas.Ab = R(:,1+Xfeas.nGc:end-1);
Xfeas.b = R(:,end);

figure;
plot(Xfeas,'b',0.1)

% cartesian product for 3D plot of explicit control law
% IOs = [ 1 zeros(1,nz-1) ];
% XoU = hybZono( [ Xfeas.Gc ; IOs*Zst.Gc ] , [ Xfeas.Gb ; IOs*Zst.Gb ] , [ Xfeas.c ; IOs*Zst.c ] , ...
% 				Xfeas.Ac , Xfeas.Ab , Xfeas.b );



%% Support function testing
% rng('default')
% Test = randomSet(2,'hybZono',2,5,3,2)'
% figure; hold on
% plot(Test,'b',0.1)
% [s,x] = supportFunc(Test,[1;0])
% plot(boundingBox(Test),'r',0.1)