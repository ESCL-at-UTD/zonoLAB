% mpt_demo_lti4
% A = [1 1; 0 1];
% A = [0.9 1; -0.5 0.9];
% B = [1; 0.5];
% C = [1 0];
% D = 0;
rng(1)
nx = 6;  % Did not work with nx = 6 nu = 2
nu = 2;
A = rand(nx,nx);
B = rand(nx,nu);
C = [1 zeros(1,nx-1)];
D = zeros(1,nu);
% A = [1 1 0; 0 1 1; 0 0 1];
% B = [1; 0.5; 0.5];
% C = [1 0 0];
% D = 0;
model = LTISystem('A', A, 'B', B, 'C', C, 'D', D);

% Next we define an MPC controller:
horizon = 5;
ctrl = MPCController(model, horizon);

% Specify the MPC problem (constraints and penalties):
% ctrl.model.x.min = [-5; -5];
% ctrl.model.x.max = [5; 5];
% ctrl.model.x.min = [-5; -5; - 5];
% ctrl.model.x.max = [5; 5; 5];
% ctrl.model.u.min = -1;
% ctrl.model.u.max = 1;
ctrl.model.x.min = -5*ones(nx,1);
ctrl.model.x.max = 5*ones(nx,1);
ctrl.model.u.min = -1*ones(nu,1);
ctrl.model.u.max = ones(nu,1);
% Use quadratic cost function
ctrl.model.x.penalty = QuadFunction(eye(model.nx));
ctrl.model.u.penalty = QuadFunction(eye(model.nu));

% Finally, we can convert the controller to an explicit form:
disp('Generating the explicit solution:');
tStart = tic;
expctrl = ctrl.toExplicit();
toc(tStart)

% % Compare optimal solutions
% x0 = [-4; 0];
% disp('Optimal control input obtained by evaluating the explicit solution:');
% Uexp = expctrl.evaluate(x0)
% 
% disp('Optimal control input obtained by solving the optimization problem on-line:');
% Uonl = ctrl.evaluate(x0)
% 
% % plot the explicit optimizer
% figure
% expctrl.feedback.fplot();
% title('PWA representation of the optimal control input')
% 
% % plot the value function
% figure
% expctrl.cost.fplot();
% title('Explicit cost function');

% plot the regions
% figure
% expctrl.partition.plot()
% title('Regions of the polyhedral partition');

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
% nh
H = feasSet.H(:,1:nz);
S = feasSet.H(:,1+nz:end-1);
f = feasSet.H(:,end);
nh = size(H,1);

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

% Yalmip was not able to find dual bounds
% z_ = sdpvar(nz,1);
% x_ = sdpvar(nx,1);
% [Constraints,details] = kkt([H*z_ + S*x_ <= f, ctrl.model.x.min <= x_ <= ctrl.model.x.max],1/2*z_'*Q*z_+q'*z_,x_);
% details.dualbounds


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

% formulate the hybrid zonotopes Z* and Xfeas
% Compute bounds
ZX = conZono([H S f]);
ZX = and(ZX,X,[zeros(nx,nz) eye(nx)]);

ZX_mapped = [H S]*ZX + zono(-f);
ZX_mappedZono = zono(ZX_mapped.Gc,ZX_mapped.c); % Conver to zonotope overapprox set and reduce computation time for boundingBox
ZX_mapped_box = boundingBox(ZX_mapped);
M = abs(ZX_mapped_box.c - ZX_mapped_box.G*ones(nh,1));

muMaxT0 = 1e2*ones(nh,1); % Initial guess smaller than it should be

% New idea: Guess large value of muMaxT, compute hybrid zonotope, call
% getLeaves. This will give us all the combinations of active constraints.
% Then use this to solve for dual varaibles as a function of parameters (as
% done in references). Compute bounds on each. If any bound is equal to
% origninal bound, increase bound and try again. Once all bounds are
% smaller than initial guess, update bounds and recompute final hybrid
% zonotope.
maxIter = 10;
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

    opts = solverOptions;
    leaves = getLeaves(Zall,opts);
    nLeaves = size(leaves,2);
    
    muMaxT = zeros(nh,1);
    for iterLeaf = 1:nLeaves
        activeIndx = find(leaves(:,iterLeaf)==1);
        if isempty(activeIndx)
            continue
        end
        H_act = H(activeIndx,:);
        S_act = S(activeIndx,:);
        f_act = f(activeIndx,:);

%         dependentIndx = find(abs(diag(1./(vecnorm(H_act').^2))*H_act*H_act')-eye(size(H_act,1))==1);
%         if isempty(dependentIndx)
%             muSet = (H_act*(Q\(H_act')))\eye(size(H_act,1))*(S_act*zono([X.Gc X.Gb],X.c)+zono(-H_act*(Q\q)-f_act)); % Does not consider X as hybZono
%             muSetBox = boundingBox(muSet);
%             ub_Box = muSetBox.c+muSetBox.G*ones(length(activeIndx),1);
%             muMaxT(activeIndx) = max(muMaxT(activeIndx),ub_Box);
%             if max(muMaxT(activeIndx)) >= 1e6
%                 disp('Full Rank. Huge bound!')
%             end
%         else 
%             if length(dependentIndx) > 2
%                 disp('More than two dependent rows!')
%             else
%             [row,col] = ind2sub([size(H_act,1) size(H_act,1)],dependentIndx);
%             end
%             for count = 1:2
%                 activeIndx0 = activeIndx;
%                 activeIndx0(row(count)) = [];
%                 H_act = H(activeIndx0,:);
%                 S_act = S(activeIndx0,:);
%                 f_act = f(activeIndx0,:);
%                 muSet = (H_act*(Q\(H_act')))\eye(size(H_act,1))*(S_act*zono([X.Gc X.Gb],X.c)+zono(-H_act*(Q\q)-f_act)); % Does not consider X as hybZono
%                 muSetBox = boundingBox(muSet);
%                 ub_Box = muSetBox.c+muSetBox.G*ones(length(activeIndx0),1);
%                 muMaxT(activeIndx0) = max(muMaxT(activeIndx0),ub_Box);
%                 if max(ub_Box) >= 1e6
%                     disp('Low Rank. Huge bound!')
%                 end
%             end
%         end
%     end

        [ Qr , Rr , E ] = qr(H_act);
        if ~isvector(Rr)
            diagr = abs(diag(Rr));
        else
            diagr = abs(Rr(1)); % Added abs
        end
        %Rank estimation
        tol = 1e-2;  % Tolerance - need to look at this more
        r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
        if r == size(H_act,1)
            dependentIndx = find(abs(diag(1./(vecnorm(H_act').^2))*H_act*H_act')-eye(size(H_act,1))==1);
            if isempty(dependentIndx)
                muSet = (H_act*(Q\(H_act')))\eye(size(H_act,1))*(S_act*zono([X.Gc X.Gb],X.c)+zono(-H_act*(Q\q)-f_act)); % Does not consider X as hybZono
                muSetBox = boundingBox(muSet);
                ub_Box = muSetBox.c+muSetBox.G*ones(length(activeIndx),1);
                muMaxT(activeIndx) = max(muMaxT(activeIndx),ub_Box);
                if max(muMaxT(activeIndx)) >= 1e6
                    disp('Full Rank. Huge bound!')
                end
            else
                if length(dependentIndx) > 2
                    disp('More than two dependent rows!')
                else
                    [row,col] = ind2sub([size(H_act,1) size(H_act,1)],dependentIndx);
                end
                for count = 1:2
                    activeIndx0 = activeIndx;
                    activeIndx0(row(count)) = [];
                    H_act = H(activeIndx0,:);
                    S_act = S(activeIndx0,:);
                    f_act = f(activeIndx0,:);
                    muSet = (H_act*(Q\(H_act')))\eye(size(H_act,1))*(S_act*zono([X.Gc X.Gb],X.c)+zono(-H_act*(Q\q)-f_act)); % Does not consider X as hybZono
                    muSetBox = boundingBox(muSet);
                    ub_Box = muSetBox.c+muSetBox.G*ones(length(activeIndx0),1);
                    muMaxT(activeIndx0) = max(muMaxT(activeIndx0),ub_Box);
                    if max(ub_Box) >= 1e6
                        disp('Low Rank. Huge bound!')
                    end
                end
            end
        else
            R1 = Rr/E;
            R1(end,:) = [];
            f12 = Qr\f_act;
            S12 = Qr\S_act;
            f2 = f12(end,:);
            S2 = S12(end,:);
            if sum(abs(S2)) <= 1e-12
                if sum(abs(f2)) <= 1e-12 
                    disp('Dont know how to handle this case.')
                else 
                    disp('Lower dimenional region not worth exploring.')
                end
            else
%                 iterLeaf
            end
            
%             indT = sort(E(1:r));
%             activeIndx = activeIndx(indT);
%             H_act = H(activeIndx,:);
%             S_act = S(activeIndx,:);
%             f_act = f(activeIndx,:);
%             muSet = (H_act*(Q\(H_act')))\eye(size(H_act,1))*(S_act*zono([X.Gc X.Gb],X.c)+zono(-H_act*(Q\q)-f_act)); % Does not consider X as hybZono
%             muSetBox = boundingBox(muSet);
%             ub_Box = muSetBox.c+muSetBox.G*ones(length(activeIndx),1);
%             muMaxT(activeIndx) = max(muMaxT(activeIndx),ub_Box);
%             if max(ub_Box) >= 1e6
%                 disp('Low Rank. Huge bound!')
%             end
        end
    end
    if min(muMaxT <= muMaxT0) == 1
        disp(['Converged in ', num2str(iter), ' iterations.'])
        break
    else
        indxIncrease = find(muMaxT >= muMaxT0);
        muMaxT0(indxIncrease) = 2*muMaxT(indxIncrease);
    end
    if iter == maxIter
        disp('Max iter reached.')
    end
end

Zmu0 = zono(diag(muMaxT)/2,muMaxT/2);
Zp = hybZono([],0.5*eye(nh),0.5*ones(nh,1),[],[],[]);

Zall = cartProd(Zb,X);
Zall = cartProd(Zall,Zmu0);
Zall = cartProd(Zall,Zp);
Zall.Ac = [Zall.Ac; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gc];
Zall.Ab = [Zall.Ab; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gb];
Zall.b = [Zall.b; -q-[Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.c];

Zall = halfspaceIntersection(Zall,[H S zeros(nh) zeros(nh)],f);
Zall = halfspaceIntersection(Zall,[-H -S zeros(nh) diag(M)],M-f);
Zall = halfspaceIntersection(Zall,[zeros(nh,nz) zeros(nh,nx) eye(nh) -diag(muMaxT)],zeros(nh,1));

% If upper-bound on mu is zero, then we can remove continue and binary
% factors cooresponding to those constraints
zeroIndx = find(muMaxT == 0);
Zall.Gc(:,nz+nx+zeroIndx) = [];
Zall.Ac(:,nz+nx+zeroIndx) = [];
Zall.c = Zall.c - Zall.Gb(:,zeroIndx)*ones(length(zeroIndx),1);
Zall.Gb(:,zeroIndx) = [];
Zall.b = Zall.b + Zall.Ab(:,zeroIndx)*ones(length(zeroIndx),1);
Zall.Ab(:,zeroIndx) = [];

Xfeas = projection(Zall,[nz+1:nz+nx]);

toc(tStart)

figure;
optPlot = plotOptions('FaceColor','b','FaceAlpha',0.1,'Display','individual');
tic
% plot(Xfeas,optPlot)
plot(projection(Xfeas,1:2),optPlot)
toc


%%
figure;
XfeasConZono = conZono([Xfeas.Gc Xfeas.Gb],Xfeas.c,[Xfeas.Ac Xfeas.Ab],Xfeas.b);
figure;
optPlot = plotOptions('FaceColor','b','FaceAlpha',1,'Display','individual');
tic
plot(XfeasConZono,optPlot)
toc

%%
figure;
ZallConZono = conZono([Zall.Gc Zall.Gb],Zall.c,[Zall.Ac Zall.Ab],Zall.b);
figure;
optPlot = plotOptions('FaceColor','b','FaceAlpha',1,'Display','individual');
tic
plot(projection(ZallConZono,nz+nx+[1:3]),optPlot)
toc


%% Next step is testing if explicit MPC and hybrid zono give the same optimal input for a range of measured states