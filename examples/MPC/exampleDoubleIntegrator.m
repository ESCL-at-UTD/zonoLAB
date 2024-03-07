% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example:
%   Representation of MPC policy for double integrator system.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Parameters
nx = 2;  % Number of states
nu = 1;  % Number of inputs
nh = 5;  % Prediction horizon
xMin = -5*ones(nx,1);  % State lower bounds
xMax =  5*ones(nx,1);  % State upper bounds
uMin = -1*ones(nu,1);  % Input lower bounds
uMax =  1*ones(nu,1);  % Input upper bounds
Q = eye(nx); % State weight matrix
R = eye(nu); % Input weight matrix
computeDualBoundsFlag = 1; % Use automated process to determine dual bounds (requires Gurobi)
includeTerminalFlag = 1;   % Use terminal penalty and constraints (requires MPT3)

% System dynamics (discrete-time double integrator)
A = [1 1; 0 1];
B = [1; 0.5];

% MPC formulation using MATLAB problem-based optimization definition

prob = optimproblem('ObjectiveSense','minimize');
u_ = optimvar('u',nu,nh);
x0_ = optimvar('x0',nx,1);

prob.Objective = 0;
prob.Constraints.Upper = optimconstr(1);
prob.Constraints.Lower = optimconstr(1);
x_ = x0_;
for k = 1:nh
    prob.Objective = prob.Objective + x_'*Q*x_ + u_(:,k)'*R*u_(:,k);
    prob.Constraints.Upper = [prob.Constraints.Upper; x_ <= xMax];
    prob.Constraints.Lower = [prob.Constraints.Lower; x_ >= xMin];
    prob.Constraints.Upper = [prob.Constraints.Upper; u_(:,k) <= uMax];
    prob.Constraints.Lower = [prob.Constraints.Lower; u_(:,k) >= uMin];
    x_ = A*x_ + B*u_(:,k);
end
prob.Constraints.Upper = [prob.Constraints.Upper; x_ <= xMax];
prob.Constraints.Lower = [prob.Constraints.Lower; x_ >= xMin];

% Optional inclusion of terminal cost and constraints
if includeTerminalFlag
    prob = addTerminal(prob,A,B,Q,R,x_,nx,xMax,xMin,uMax,uMin); % Local function below
end

% Identify QP matrices
problem = prob2struct(prob);

% Current formulation only allows for inequality constraints. 
% Reformulating equality constraints as inequality constraints.
A_con = full([problem.Aineq; problem.Aeq; -problem.Aeq]);
B_con = full([problem.bineq; problem.beq; -problem.beq]);

% Specify strcture for multi-parametric Quadratic Program (mpQP)
% with objective function and constraints in the form:
% min_U 1/2U'*Q*U + x'*P*U + q'*U
%	s.t. HU + Sx <= f
mpQP.nz = nu*nh; % Number of decision variables
mpQP.Q = full(problem.H(1:mpQP.nz,1:mpQP.nz));
mpQP.P = full(problem.H(1+mpQP.nz:end,1:mpQP.nz));
mpQP.q = full(problem.f(1:mpQP.nz)')';
mpQP.H = full(A_con(:,1:mpQP.nz));
mpQP.S = full(A_con(:,mpQP.nz+1:end));
mpQP.f = full(B_con);
mpQP.nx = size(mpQP.S,2);  % Number of parameters
mpQP.nh = size(mpQP.S,1);  % Number of inequality constraints

% Define bounded set for parameters 
mpQP.X = zono(diag((xMax-xMin)/2),(xMax+xMin)/2);

% Compute or specify upper bounds on dual variables (muMax)
if computeDualBoundsFlag
    muMax = [];
else
    muMax = 1e4*ones(mpQP.nh,1);
end

% Compute optimal control law map
tStart = tic;
Zstar = mpQPMap(mpQP,muMax);
tMPC = toc(tStart)

% Plot control law
optPlot = plotOptions('FaceColor','b','FaceAlpha',0.1,'Display','individual');
figure;
tStart = tic;
plot(Zstar,optPlot) % 27 leaves are plotted but 2 are only line segments
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
zlabel('$u^*$','interpreter','latex')
view(-130,35)
set(gca,'fontsize',18,'fontname','times new roman')
tPlot = toc(tStart)

function prob = addTerminal(prob,A,B,Q,R,x_,nx,xMax,xMin,uMax,uMin)
    % Compute terminal cost
    [K, P] = dlqr(A,B,Q,R);
    K = -K; % u = K*x
    prob.Objective = prob.Objective + x_'*P*x_;
    % Compute terminal constraint as maximumal positive invariant set
    AK = A+B*K; % Closed-loop LQR System: x^+ = (A+B*K)x
    Omega.H = [eye(nx); -eye(nx); K; -K];
    Omega.f = [xMax;-xMin;uMax;-uMin];
    OmegaMPTOld = Polyhedron('A',Omega.H,'b',Omega.f);
    for indx = 1:1e3 % Iteration limit
        Omega.H = [Omega.H*AK;Omega.H];
        Omega.f = [Omega.f;Omega.f];
        OmegaMPT = Polyhedron('A',Omega.H,'b',Omega.f);
        OmegaMPT.minHRep;
        if OmegaMPT == OmegaMPTOld
            break
        else
            OmegaMPTOld = OmegaMPT;
            Omega.H = OmegaMPT.H(:,1:end-1);
            Omega.f = OmegaMPT.H(:,end);
        end
    end
    prob.Constraints.Terminal = optimconstr(1);
    prob.Constraints.Terminal = Omega.H*x_ <= Omega.f;
end