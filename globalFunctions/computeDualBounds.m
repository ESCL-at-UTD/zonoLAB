% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Compute upper bounds on the dual variables for a multi-parametric 
%       Quadratic Program (mpQP)
%   Syntax:
%       muMax = computeDualBounds(mpQP,Zb,M)
%   Inputs:
%       mpQP - structure containing all mpQP data/matrices
%       Zb   - zono object bounding mpQP decision variable
%       M    - nh x 1 vector bound such that -M <= Hz + Sx - f for all z and x
%   Outputs:
%       muMax - nh x 1 vector bound such that 0 <= mu <= muMax
%   Notes:
%       This code is a work in progress.
%       nh is the number of inequality constrains in the mpQP.
%       General idea: Guess value of muMax, compute hybrid zonotope, call
%       getLeaves. This gives all the combinations of active constraints.
%       Then use to solve for dual varaibles as a function of parameters (as
%       done in Bemporad et. al. 2002). Compute bounds on each. If any bound 
%       is equal toorigninal bound, increase bound and try again. Once all 
%       bounds are smaller than initial guess, update bounds and recompute 
%       final hybrid zonotope.
%       Bemporad, Morari, Dua, and Pistikopoulos, "The explicit linear
%       quadratic regulator for constrained systems," Automatica, 2002.
%       Section 4.1.1. Degeneracy
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function muMax = computeDualBounds(mpQP,Zb,M)

% Unpack problem data
nz = mpQP.nz;
nx = mpQP.nx;
nh = mpQP.nh;
Q  = mpQP.Q;
P  = mpQP.P;
q  = mpQP.q;
H  = mpQP.H;
S  = mpQP.S;
f  = mpQP.f;
X  = mpQP.X;

muMax0 = 1e2*ones(nh,1); % Initial guess smaller than it should be

maxIter = 10;            % Maximum iterations 
for iter = 1:maxIter
    Zmu0 = zono(diag(muMax0)/2,muMax0/2);
    Zp = hybZono([],0.5*eye(nh),0.5*ones(nh,1),[],[],[]);

    Zall = cartProd(Zb,X);
    Zall = cartProd(Zall,Zmu0);
    Zall = cartProd(Zall,Zp);
    Zall.Ac = [Zall.Ac; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gc];
    Zall.Ab = [Zall.Ab; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gb];
    Zall.b = [Zall.b; -q-[Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.c];

    Zall = halfspaceIntersection(Zall,[H S zeros(nh) zeros(nh)],f);
    Zall = halfspaceIntersection(Zall,[-H -S zeros(nh) diag(M)],M-f);
    Zall = halfspaceIntersection(Zall,[zeros(nh,nz) zeros(nh,nx) eye(nh) -diag(muMax0)],zeros(nh,1));

    opts = solverOptions;
    leaves = getLeaves(Zall,opts);
    nLeaves = size(leaves,2);

    muMax = zeros(nh,1);
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
        % Rank estimation
        tol = 1e-2;        % Tolerance - need to look at this more
        r = find(diagr >= tol*diagr(1), 1, 'last'); %rank estimation
        if r == size(H_act,1)
            H_norm = diag(1./vecnorm(H_act'))*H_act;
            dependentIndx = find(abs(H_norm*H_norm')-eye(size(H_act,1))==1);
            if isempty(dependentIndx)
                muSet = (H_act*(Q\(H_act')))\eye(size(H_act,1))*(S_act*X+zono(-H_act*(Q\q)-f_act)); % Does not consider X as hybZono
                muSetBox = boundingBox(muSet);
                ub_Box = muSetBox.c+muSetBox.G*ones(length(activeIndx),1);
                muMax(activeIndx) = max(muMax(activeIndx),ub_Box);
                if max(muMax(activeIndx)) >= 1e6
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
                    muSet = (H_act*(Q\(H_act')))\eye(size(H_act,1))*(S_act*X+zono(-H_act*(Q\q)-f_act)); % Does not consider X as hybZono
                    muSetBox = boundingBox(muSet);
                    ub_Box = muSetBox.c+muSetBox.G*ones(length(activeIndx0),1);
                    muMax(activeIndx0) = max(muMax(activeIndx0),ub_Box);
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
                muSet = (R1*(Q\(R1')))\eye(r)*(S1*X+zono(-R1*(Q\q)-f1)); % Does not consider X as hybZono
                muSet = inv(Qr')*cartProd(muSet,zono(zeros(size(H_act,1)-r,1)));
                muSetBox = boundingBox(muSet);
                ub_Box = muSetBox.c+muSetBox.G*ones(length(activeIndx),1);
                muMax(activeIndx) = max(muMax(activeIndx),ub_Box);
                if max(ub_Box) >= 1e6
                    disp('Low Rank. Huge bound!')
                end
            end
        end
    end
    if min(muMax <= muMax0) == 1
        disp(['Converged in ', num2str(iter), ' iterations.'])
        break
    else
        indxIncrease = find(muMax >= muMax0);
        muMax0(indxIncrease) = 2*muMax(indxIncrease);
    end
    if iter == maxIter
        disp('Max iter reached.')
    end
end

end