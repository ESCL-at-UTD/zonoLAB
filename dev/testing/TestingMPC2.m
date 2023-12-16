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

% Find and remove all-zero rows of H (corresponding to constraints on
% parameters, which are not needed)
zeroIndx = find(sum(abs(H),2)==0);
H(zeroIndx,:) = [];
S(zeroIndx,:) = [];
f(zeroIndx,:) = [];

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
% muMaxT = 1e3*ones(nh,1);

% Compute bounds
ZX = conZono([H S f]);
ZX = and(ZX,X,[zeros(nx,nz) eye(nx)]);

ZX_mapped = [H S]*ZX + zono(-f);
ZX_mappedZono = zono(ZX_mapped.Gc,ZX_mapped.c); % Conver to zonotope overapprox set and reduce computation time for boundingBox
ZX_mapped_box = boundingBox(ZX_mapped);
M = abs(ZX_mapped_box.c - ZX_mapped_box.G*ones(nh,1));

% muMaxT0 = 1e-1*ones(nh,1); % Initial guess smaller than it should be
muMaxT0 = 1e2*ones(nh,1); % Initial guess smaller than it should be
% muMaxT0(1) = 300;
% muMaxT0(6) = 300;
% muMaxT0(13) = 1e3;
% muMaxT0(18) = 1e3;

% maxIter = 10;
maxIter = 1;
for iter = 1:maxIter
    Zmu0 = zono(diag(muMaxT0)/2,muMaxT0/2);
    Zp = hybZono([],0.5*eye(nh),0.5*ones(nh,1),[],[],[]);

    Zall = cartProd(Zb,X);
    Zall = cartProd(Zall,Zmu0);
    Zall = cartProd(Zall,Zp);
    Zall.Ac = [Zall.Ac; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gc];
    Zall.Ab = [Zall.Ab; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gb];
    Zall.b = [Zall.b; -q-[Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.c];

    Zall = halfspaceIntersection(Zall,[H S zeros(nh) zeros(nh)],f);
    Zall = halfspaceIntersection(Zall,[-H -S zeros(nh) diag(M)],M-f);
    Zall = halfspaceIntersection(Zall,[zeros(nh,nz) zeros(nh,nx) eye(nh) -diag(muMaxT0)],zeros(nh,1));

    Zmu = projection(Zall,nz+nx+[1:nh]);
    Zmu_box = boundingBox(Zmu); % This operation is slow (Maybe something
%     using getLeaves would be better
% 
%     % Find if upper-bounds on mu equal the initial value
    muMaxT = Zmu_box.c+Zmu_box.G*ones(nh,1);
%     indxMu = muMaxT == muMaxT0;
% %     indxMu
%     if min(muMaxT) == 0 % Likely sign that bounds were too small
%         muMaxT = 10*muMaxT0;
%     else
% %         muMaxT(indxMu) = 2*muMaxT0(indxMu);
%     end
% 
% %     for i = 1:nh
% %         if muMaxT(i) == 0
% %             muMaxT(i) = 2*muMaxT0(i); % Increase bound
% %         else
% %             if indxMu(i) == 1
% %                 if min(H(i,:) == zeros(1,nz))  % Not a problem if corresponding row of H is all zeros
% %                     indxMu(i) = 0;
% %                     muMaxT(i) = 0;
% %                 else
% %                     muMaxT(i) = 2*muMaxT0(i); % Increase bound
% %                 end
% %             else
% %                 muMaxT(i) = Zmu_box.c(i)+Zmu_box.G(i,:)*ones(nh,1); % Set bound to computed value
% %             end
% %         end
% %     end
% %     iter
%     if min(muMaxT == muMaxT0) == 1
%         disp('Converged.')
%         break
%     end
% %     muMaxT
%     muMaxT0 = muMaxT;
end


% New idea: Guess large value of muMaxT, compute hybrid zonotope, call
% getLeaves. This will give us all the combinations of active constraints.
% Then use this to solve for dual varaibles as a function of parameters (as
% done in references). Compute bounds on each. If any bound is equal to
% origninal bound, increase bound and try again. Once all bounds are
% smaller than initial guess, update bounds and recompute final hybrid
% zonotope.


% H_aligned = zeros(nh);
% for i = 1:nh
%     for j = 1:i-1
%         H_aligned(i,j) = H(i,:)*H(j,:)'/(norm(H(i,:)*norm(H(j,:))));
%     end
% end
% H_aligned = H_aligned == 1;

% Find if upper-bounds on mu equal the initial value
% indxMu = Zmu_box.c+Zmu_box.G*ones(nh,1) == muMaxT0;
% 
% muMaxT = Zmu_box.c+Zmu_box.G*ones(nh,1);
% 
% for i = 1:nh
%     if indxMu(i) == 1
%         H(i,:)
%         if min(H(i,:) == zeros(1,nz))  % Not a problem if corresponding row of H is all zeros
%             indxMu(i) = 0;
%             muMaxT(i) = 0;
%         end
%     end
% end
% if sum(indxMu) == 0
%     disp('All mu are bounded properly.')
% else
%     error('Some mu are not bounded properly.')
% end
% 
% 
% Zmu = zono(diag(muMaxT)/2,muMaxT/2);
% 
% Zall = cartProd(Zb,X);
% Zall = cartProd(Zall,Zmu);
% Zall = cartProd(Zall,Zp);
% Zall.Ac = [Zall.Ac; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gc];
% Zall.Ab = [Zall.Ab; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gb];
% Zall.b = [Zall.b; -q-[Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.c];
% 
% Zall = halfspaceIntersection(Zall,[H S zeros(nh) zeros(nh)],f);
% Zall = halfspaceIntersection(Zall,[-H -S zeros(nh) diag(M)],M-f);
% Zall = halfspaceIntersection(Zall,[zeros(nh,nz) zeros(nh,nx) eye(nh) -diag(muMaxT0)],zeros(nh,1));

Xfeas = projection(Zall,[nz+1:nz+nx]);

% R = rref([Xfeas.Ac Xfeas.Ab Xfeas.b]);
% Xfeas.Ac = R(:,1:Xfeas.nGc);
% Xfeas.Ab = R(:,1+Xfeas.nGc:end-1);
% Xfeas.b = R(:,end);

toc(tStart)

figure;
plot(Xfeas,'b',0.1)
