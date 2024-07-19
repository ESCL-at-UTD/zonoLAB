% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Global function:
%       Returns a hybrid/constrained zonotope based on a collection of nV vertices
%   Syntax:
%       Z = vPoly2Zono(V,M)
%   Inputs:
%       V - n x nV matrix with each column defining a vertex in R^n
%       M - (optional) nV x N matrix indicating which of the N sets the nV vertices belong
%           if M is not given, all vertices belong to the same contiguous set and a conZono is returned
%   Outputs:
%       Z - hybrid or constrained zonotope in R^n
%   Notes:
%       Hybrid/constrained zonotope represents a set based on the nV vertices stored
%       in V and the incidence matrix M
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = vPoly2Zono(V,varargin)
    [n,nV] = size(V);

    if nargin == 1
        N = 1;
        M = ones(nV,1);
    else
        M = varargin{1};
        N = size(M,2);
    end
    
    Gc = 0.5*[eye(nV);zeros(N,nV)];
    Gb = 0.5*[zeros(nV,N);eye(N)];
    c = 0.5*ones(nV+N,1);
    Ac = [ones(1,nV);zeros(1,nV)];
    Ab = [zeros(1,N);ones(1,N)];
    b = [2-nV;2-N];
    
    Q = hybZono(Gc,Gb,c,Ac,Ab,b);
    
    R = [eye(nV) -M];
    D = halfspaceIntersection(Q,eye(nV),zeros(nV,1),R);
    
    Z = [V zeros(n,N)]*D;

    if N == 1
        leaf = getLeaves(Z,{});
        out = conZono(Z.Gc,Z.c+Z.Gb*leaf,Z.Ac,Z.b-Z.Ab*leaf);
    else
        out = Z;
    end