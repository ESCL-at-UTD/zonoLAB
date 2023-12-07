% mpt_demo_lti4
% A = [1 1; 0 1];
A = [0.9 1; -0.5 0.9];
B = [1; 0.5];
nx = 2;  % Did not work with nx = 6 nu = 2 (seems to be working now)
nu = 1;
% rng(1)
% A = rand(nx,nx);
% B = rand(nx,nu);
C = [1 zeros(1,nx-1)];
D = zeros(1,nu);
model = LTISystem('A', A, 'B', B, 'C', C, 'D', D);

% Next we define an MPC controller:
horizon = 10;
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

% % % Compare optimal solutions
% x0 = [-5; -1];
% disp('Optimal control input obtained by evaluating the explicit solution:');
% Uexp = expctrl.evaluate(x0)
% 
% disp('Optimal control input obtained by solving the optimization problem on-line:');
% Uonl = ctrl.evaluate(x0)
% 
% plot the explicit optimizer
figure
expctrl.feedback.fplot();
title('PWA representation of the optimal control input')
% % 
% % % plot the value function
% % figure
% % expctrl.cost.fplot();
% % title('Explicit cost function');
% 
% % plot the regions
% figure
% expctrl.partition.plot()
% title('Regions of the polyhedral partition');

%% Matlab problem-based optimization problem definition
% prob = optimproblem('ObjectiveSense','minimize');
% x_ = optimvar('x',nx,horizon+1,'LowerBound',repmat(ctrl.model.x.min,1,horizon+1),'UpperBound',repmat(ctrl.model.x.max,1,horizon+1));
% u_ = optimvar('u',nu,horizon,'LowerBound',repmat(ctrl.model.u.min,1,horizon),'UpperBound',repmat(ctrl.model.u.max,1,horizon));
% 
% prob.Objective = 0;
% for k = 1:horizon
%     prob.Objective = prob.Objective + x_(:,k)'*eye(nx)*x_(:,k) + u_(:,k)'*eye(nu)*u_(:,k);
%     prob.Constraints.dynamics = x_(:,k+1) == A*x_(:,k) + B*u_(:,k);
% end
% x0.x = [-5;-1];
% problem = prob2struct(prob,[-5;-1]);
% 
% 
% u_star = quadprog(problem)

prob = optimproblem('ObjectiveSense','minimize');
u_ = optimvar('u',nu,horizon);
x0_ = optimvar('x0',nx,1);


prob.Objective = 0;
prob.Constraints.Upper = optimconstr(1);
prob.Constraints.Lower = optimconstr(1);
x_ = x0_;
for k = 1:horizon
    prob.Objective = prob.Objective + x_'*eye(nx)*x_ + u_(:,k)'*eye(nu)*u_(:,k);
    prob.Constraints.Upper = [prob.Constraints.Upper; x_ <= ctrl.model.x.max];
    prob.Constraints.Lower = [prob.Constraints.Lower; x_ >= ctrl.model.x.min];
    prob.Constraints.Upper = [prob.Constraints.Upper; u_(:,k) <= ctrl.model.u.max];
    prob.Constraints.Lower = [prob.Constraints.Lower; u_(:,k) >= ctrl.model.u.min];
    x_ = A*x_ + B*u_(:,k);
end
prob.Constraints.Upper = [prob.Constraints.Upper; x_ <= ctrl.model.x.max];
prob.Constraints.Lower = [prob.Constraints.Lower; x_ >= ctrl.model.x.min];

problem = prob2struct(prob);

A_con = full([problem.Aineq; problem.Aeq; -problem.Aeq]);
B_con = full([problem.bineq; problem.beq; -problem.beq]);

nz = nu*horizon;

Matrices.G = A_con(:,1:nz);
Matrices.E = -A_con(:,nz+1:end);
Matrices.W = B_con;
Matrices.H = full(problem.H(1:nz,1:nz));
Matrices.F = full(problem.H(1+nz:end,1:nz));
Matrices.Cf = full(problem.f(1:nz)');

%% zonoLab Implementation
%   min_U 1/2 U'*H*U + x'*F*U + Cf*U + x'*Y*x + Cx*x + Cc
%    s.t. G*U <= W + E*x

% Matrices = getMatrices(ctrl);

% dimensions of the mp-QP
% nz = size(Matrices.H,2);
% nx = size(Matrices.E,2);

% represent the mp-QP as a hybrid zonotope

tStart = tic;

% define parameter set over all of the bounded parameters
Gx = diag((ctrl.model.x.max-ctrl.model.x.min)/2);
cx = (ctrl.model.x.max+ctrl.model.x.min)/2;
X = zono(Gx,cx);
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
H = feasSet.H(:,1:nz);
S = feasSet.H(:,1+nz:end-1);
f = feasSet.H(:,end);
nh = size(H,1);

% H_norm = diag(1./vecnorm(H'))*H;
% dependentIndx = find(abs(H_norm*H_norm')-eye(nh)==1);
% [row,col] = ind2sub([nh nh],dependentIndx);
% pairs = [row,col];
% indRemove = find(pairs(:,1)-pairs(2) <=0 );
% pairs(indRemove,:)= [];
% nPairs = size(pairs,1);
% pairsMatrix = zeros(nPairs,nh);
% pairsMatrix(sub2ind([nPairs,nh],[1:nPairs]',pairs(:,1))) = 1;
% pairsMatrix(sub2ind([nPairs,nh],[1:nPairs]',pairs(:,2))) = 1;

% now convert to standard form
% min_z 1/2z'*Q*z + q'*z
%	s.t. Hz + Sx <= f
%	where z = U + inv(Q)*P'*x 
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

% formulate the hybrid zonotopes Z* and Xfeas
% Compute bounds
ZX = conZono([H S f]);
ZX = and(ZX,X,[zeros(nx,nz) eye(nx)]);

ZX_mapped = [H S]*ZX + zono(-f);
ZX_mappedZono = zono(ZX_mapped.Gc,ZX_mapped.c); % Conver to zonotope overapprox set and reduce computation time for boundingBox
ZX_mapped_box = boundingBox(ZX_mappedZono);
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
%     Zall = halfspaceIntersection(Zall,[zeros(1,nz) zeros(1,nx) zeros(1,nh) ones(1,nh)],nz);
%     Zall = halfspaceIntersection(Zall,[zeros(nPairs,nz) zeros(nPairs,nx) zeros(nPairs,nh) pairsMatrix],ones(nPairs,1)); % Does not seem to help speed things up

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
            H_norm = diag(1./vecnorm(H_act'))*H_act;
            dependentIndx = find(abs(H_norm*H_norm')-eye(size(H_act,1))==1);
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
            R12 = Rr/E;
            R1 = R12(1:r,:);
            R2 = R12(r+1:end,:);
            f12 = Qr\f_act;
            S12 = Qr\S_act;
            f1 = f12(1:r,:);
            S1 = S12(1:r,:);
            f2 = f12(r+1:end,:);
            S2 = S12(r+1:end,:);
            if sum(sum(abs(S2))) <= 1e-12
                if sum(sum(abs(f2))) <= 1e-12 
                    disp('Dont know how to handle this case.')
                else 
                    disp('Lower dimenional region not worth exploring.')
                end
            else
%                 iterLeaf
                muSet = (R1*(Q\(R1')))\eye(r)*(S1*zono([X.Gc X.Gb],X.c)+zono(-R1*(Q\q)-f1)); % Does not consider X as hybZono
                muSet = inv(Qr')*cartProd(muSet,zono(zeros(size(H_act,1)-r,1)));
                muSetBox = boundingBox(muSet);
                ub_Box = muSetBox.c+muSetBox.G*ones(length(activeIndx),1);
                muMaxT(activeIndx) = max(muMaxT(activeIndx),ub_Box);
                if max(ub_Box) >= 1e6
                    disp('Low Rank. Huge bound!')
                end
            end
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

% Zmu0 = zono(diag(muMaxT)/2,muMaxT/2);
% Zp = hybZono([],0.5*eye(nh),0.5*ones(nh,1),[],[],[]);
% 
% Zall = cartProd(Zb,X);
% Zall = cartProd(Zall,Zmu0);
% Zall = cartProd(Zall,Zp);
% Zall.Ac = [Zall.Ac; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gc];
% Zall.Ab = [Zall.Ab; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gb];
% Zall.b = [Zall.b; -q-[Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.c];
% 
% Zall = halfspaceIntersection(Zall,[H S zeros(nh) zeros(nh)],f);
% Zall = halfspaceIntersection(Zall,[-H -S zeros(nh) diag(M)],M-f);
% Zall = halfspaceIntersection(Zall,[zeros(nh,nz) zeros(nh,nx) eye(nh) -diag(muMaxT)],zeros(nh,1));
% % Zall = halfspaceIntersection(Zall,[zeros(1,nz) zeros(1,nx) zeros(1,nh) ones(1,nh)],nz);

% If upper-bound on mu is zero, then we can remove continuous and binary
% factors cooresponding to those constraints
% zeroIndx = find(muMaxT == 0);
% Zall.Gc(:,nz+nx+zeroIndx) = [];
% Zall.Ac(:,nz+nx+zeroIndx) = [];
% Zall.c = Zall.c - Zall.Gb(:,zeroIndx)*ones(length(zeroIndx),1);
% Zall.Gb(:,zeroIndx) = [];
% Zall.b = Zall.b + Zall.Ab(:,zeroIndx)*ones(length(zeroIndx),1);
% Zall.Ab(:,zeroIndx) = [];
% Look at this again. Should be able to remove some rows too.  May be best
% to identify zeroIndx first, redine matrices and then construct Zall.

zeroIndx = find(muMaxT == 0);
nh = nh - length(zeroIndx);
muMaxTRed = muMaxT;
muMaxTRed(zeroIndx) = [];
HRed = H;
HRed(zeroIndx,:) = [];
SRed = S;
SRed(zeroIndx,:) = [];
fRed = f;
fRed(zeroIndx) = [];
MRed = M;
MRed(zeroIndx) = [];

Zmu0 = zono(diag(muMaxTRed)/2,muMaxTRed/2);
Zp = hybZono([],0.5*eye(nh),0.5*ones(nh,1),[],[],[]);

Zall = cartProd(Zb,X);
Zall = cartProd(Zall,Zmu0);
Zall = cartProd(Zall,Zp);
Zall.Ac = [Zall.Ac; [Q zeros(nz,nx) HRed' zeros(nz,nh)]*Zall.Gc];
Zall.Ab = [Zall.Ab; [Q zeros(nz,nx) HRed' zeros(nz,nh)]*Zall.Gb];
Zall.b = [Zall.b; -q-[Q zeros(nz,nx) HRed' zeros(nz,nh)]*Zall.c];

Zall = halfspaceIntersection(Zall,[HRed SRed zeros(nh) zeros(nh)],fRed);
Zall = halfspaceIntersection(Zall,[-HRed -SRed zeros(nh) diag(MRed)],MRed-fRed);
Zall = halfspaceIntersection(Zall,[zeros(nh,nz) zeros(nh,nx) eye(nh) -diag(muMaxTRed)],zeros(nh,1));
Zall = halfspaceIntersection(Zall,[zeros(1,nz) zeros(1,nx) zeros(1,nh) ones(1,nh)],nz);

Xfeas = projection(Zall,[nz+1:nz+nx]);

toc(tStart)

optPlot = plotOptions('FaceColor','b','FaceAlpha',0.1,'Display','individual');
% figure;
% tic
% % plot(Xfeas,optPlot)
% plot(projection(Xfeas,1:2),optPlot)
% drawnow
% toc

% Control input U = z - inv(Q)*P'*x
ZallUall = [eye(nz+nx+2*nh); [eye(nz) -inv(Q)*P' zeros(nz,nh) zeros(nz,nh)]]*Zall;
% Zstar = projection(ZallUall,[nz+1 nz+2 nz+nx+2*nh+1]);
Zstar = projection(ZallUall,[nz+[1:nx] nz+nx+2*nh+1]);
figure;
tic
plot(Zstar,optPlot)
toc


%%
% figure;
% XfeasConZono = conZono([Xfeas.Gc Xfeas.Gb],Xfeas.c,[Xfeas.Ac Xfeas.Ab],Xfeas.b);
% figure;
% optPlot = plotOptions('FaceColor','b','FaceAlpha',1,'Display','individual');
% tic
% plot(XfeasConZono,optPlot)
% toc

%%
% figure;
% ZallConZono = conZono([Zall.Gc Zall.Gb],Zall.c,[Zall.Ac Zall.Ab],Zall.b);
% figure;
% optPlot = plotOptions('FaceColor','b','FaceAlpha',1,'Display','individual');
% tic
% plot(projection(ZallConZono,nz+nx+[1:3]),optPlot)
% toc


%% Next step is testing if explicit MPC and hybrid zono give the same optimal input for a range of measured states
% Could also modify plotting code to only plot full dimensional sets to
% check if the number of sets matches MPT3

% H*1 +(S+H*inv(Q)*P')*[5;1]
% H*1 + S*[-5;-1]
% f
% 
% full(Matrices.G*1)
% Matrices.W + Matrices.E*[-5;-1]

% 
% [A,b] = equationsToMatrix(constraints,[decisions(:) ; parameters(:)]);
% % check if constraints have other terms
% if ~isempty(symvar(b))
% 	error('The constraints contain variables that are not decisions or parameters')
% end
% A = double(A);
% b = double(b);
% H = A(:,1:nz);
% Huser = H;
% S = A(:,nz+1:end);
% Suser = S;
% f = b;
% fuser = f;

%%

% Ydata = ctrl.toYALMIP();
% 
% [A,b] = equationsToMatrix(Ydata.constraints,[Ydata.variables.u; Ydata.variables.x(:,1)]);
% % check if constraints have other terms
% if ~isempty(symvar(b))
% 	error('The constraints contain variables that are not decisions or parameters')
% end
% A = double(A);
% b = double(b);
% H = A(:,1:nz);
% Huser = H;
% S = A(:,nz+1:end);
% Suser = S;
% f = b;
% fuser = f;

%%
% x = [2;-5];
% ZStar = projection(Zall,[1:nz]);
% ZStar.Ac = [ZStar.Ac; [zeros(nx,nz) eye(nx) zeros(nx,2*nh)]*Zall.Gc];
% ZStar.Ab = [ZStar.Ab; [zeros(nx,nz) eye(nx) zeros(nx,2*nh)]*Zall.Gb];
% ZStar.b = [ZStar.b; x-[zeros(nx,nz) eye(nx) zeros(nx,2*nh)]*Zall.c];
% leavesZStar = getLeaves(ZStar,opts);
% nLeavesZStar = size(leavesZStar,2);
% for i = 1:nLeavesZStar
%     Zi = conZono(ZStar.Gc,ZStar.c+ZStar.Gb*leavesZStar(:,i),ZStar.Ac,ZStar.b-ZStar.Ab*leavesZStar(:,i));
%     
%     [vi,fi] = plot(Zi,optPlot)
% end
% figure;
% [v,f] = plot(ZStar,optPlot);


