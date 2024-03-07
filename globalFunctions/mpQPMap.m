% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Compute optimal input-output map for for a multi-parametric 
%       Quadratic Program (mpQP) as a hybrid zonotope
%   Syntax:
%       Zstar = mpQPMap(mpQP,muMax)
%   Inputs:
%       mpQP  - structure containing all mpQP data/matrices
%       muMax - nh x 1 vector bound on dual variables such that 0 <= mu <= muMax
%   Outputs:
%       Zstar - hybZono object
%   Notes:
%       muMax is an optional input: can specify bounds if known or set 
%       muMax = [] and this code with compute bounds (can be slow)    
%       Zstar is a continuous piecewise-linear set that contains [x' z']'
%       where x are the mpQP parameters (inputs) and z are the
%       cooresponding optimal decision variables (outputs)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function Zstar = mpQPMap(mpQP,muMax)

% check if the problem is strictly convex
[~,p1] = chol(mpQP.Q);
if p1
	error('The objective function is not positive definite.')
end

% Find and remove all-zero rows of H (corresponding to constraints on
% parameters, which are not needed)
zeroIndx = find(sum(abs(mpQP.H),2)==0);
mpQP.H(zeroIndx,:) = [];
mpQP.S(zeroIndx,:) = [];
mpQP.f(zeroIndx,:) = [];
mpQP.nh = size(mpQP.S,1);

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

% Store Q^-1
Qinv = inv(Q);

% Convert to standard form
% min_z 1/2z'*Q*z + q'*z
%	s.t. Hz + Sx <= f
%	where z = U - inv(Q)*P'*x
S = S - H*Qinv*P';

% Find over approximative zonotope Zb > Z*
% s.t. Hz + Sx <= f
% x in X
f_max = f - S*X.c + sum(abs(S*X.G),2);
Z = conZono([H f_max]);
Zb = boundingBox(Z);

% Compute bounds such that -M <= Hz + Sx - f for all z and x
ZX = conZono([H S f]);
ZX = and(ZX,X,[zeros(nx,nz) eye(nx)]);
ZX_mapped = [H S]*ZX + zono(-f);
% Convert to zonotope overapprox set and reduce computation time for boundingBox
ZX_mappedZono = zono(ZX_mapped.G,ZX_mapped.c); 
ZX_mapped_box = boundingBox(ZX_mappedZono);
M = abs(ZX_mapped_box.c - ZX_mapped_box.G*ones(nh,1));

% Compute dual bounds such that 0 <= mu <= muMax
if ~isempty(muMax) % Dual bounds provided
    muMax(zeroIndx,:) = [];  % If any constraints were removed
else                     % Dual bounds need to be computed (which is slow)
    muMax = computeDualBounds(mpQP,Zb,M);
end

% Find if any muMax are zero (constraint is never active)
zeroMuMaxIndx = find(muMax == 0);
nh = nh - length(zeroMuMaxIndx);
muMax(zeroMuMaxIndx) = [];
H(zeroMuMaxIndx,:) = [];
S(zeroMuMaxIndx,:) = [];
f(zeroMuMaxIndx) = [];
M(zeroMuMaxIndx) = [];

Zmu0 = zono(diag(muMax)/2,muMax/2);
Zp = hybZono([],0.5*eye(nh),0.5*ones(nh,1),[],[],[]);

Zall = cartProd(Zb,X);
Zall = cartProd(Zall,Zmu0);
Zall = cartProd(Zall,Zp);
Zall.Ac = [Zall.Ac; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gc];
Zall.Ab = [Zall.Ab; [Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.Gb];
Zall.b = [Zall.b; -q-[Q zeros(nz,nx) H' zeros(nz,nh)]*Zall.c];

Zall = halfspaceIntersection(Zall,[H S zeros(nh) zeros(nh)],f);
Zall = halfspaceIntersection(Zall,[-H -S zeros(nh) diag(M)],M-f);
Zall = halfspaceIntersection(Zall,[zeros(nh,nz) zeros(nh,nx) eye(nh) -diag(muMax)],zeros(nh,1));
Zall = halfspaceIntersection(Zall,[zeros(1,nz) zeros(1,nx) zeros(1,nh) ones(1,nh)],nz);

ZallUall = [eye(nz+nx+2*nh); [eye(nz) -Qinv*P' zeros(nz,nh) zeros(nz,nh)]]*Zall;
Zstar = projection(ZallUall,[nz+[1:nx] nz+nx+2*nh+1]);

end