% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Global function:
%       Returns a hybrid zonotope based on a collection of nV vertices
%   Syntax:
%       Z = vPoly2hybZono({V,M})
%   Inputs:
%       V - n x nV matrix with each column defining a vertex in R^n
%       M - nV x N matrix indicating which of the N sets the nV vertices belong
%   Outputs:
%       Z - hybrid zonotope in R^n
%   Notes:
%       Hybrid zonotope represents a set based on the nV vertices stored
%       in V and the incidence matrix M
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = vPoly2hybZono(in)

V = in{1};
M = in{2};

n = size(V,1);
nV = size(V,2);
N = size(M,2);

Gc = 0.5*[eye(nV);zeros(N,nV)];
Gb = 0.5*[zeros(nV,N);eye(N)];
c = 0.5*ones(nV+N,1);
Ac = [ones(1,nV);zeros(1,nV)];
Ab = [zeros(1,N);ones(1,N)];
b = [2-nV;2-N];

Q = hybZono(Gc,Gb,c,Ac,Ab,b);

R = [eye(nV) -M];
D = halfspaceIntersection(Q,eye(nV),zeros(nV,1),R);

out = [V zeros(n,N)]*D;

