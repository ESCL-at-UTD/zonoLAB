% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Global function:
%       Returns a random zonotopic sets of the desired type and dimension
%   Syntax:
%       Z = randomSet(seed,type,n,nGc,nGb,nC)
%   Inputs:
%       seed - seed for random number generator
%       type - 'zono', 'conZono', or 'hybZono'
%       n    - set dimension
%       nGc  - number of continuous generators
%       nGb  - number of binary generators
%       nC   - number of constraints
%   Outputs:
%       Z - zonotopic set in R^n (hybZono, conZono, or zono object)
%   Notes:
%       Use nGc = nG and nGb = 0 when defining a zono or conZono
%       Since rand generates values from a uniform distribution in the
%       interval (0,1), 1-2*rand is used to generate value from a uniform
%       distribution in the interval (-1,1)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = randomSet(seed,type,n,nGc,nGb,nC)

rng(seed)

switch type
    case 'zono'
        G = 1-2*rand(n,nGc);
        c = 1-2*rand(n,1);
        out = zono(G,c);
    case 'conZono'
        G = 1-2*rand(n,nGc);
        c = 1-2*rand(n,1);
        A = 1-2*rand(nC,nGc);
        b = 1-2*rand(nC,1);
        out = conZono(G,c,A,b);
    case 'hybZono'
        Gc = 1-2*rand(n,nGc);
        Gb = 1-2*rand(n,nGb);
        c = 1-2*rand(n,1);
        Ac = 1-2*rand(nC,nGc);
        Ab = 1-2*rand(nC,nGb);
        b = 1-2*rand(nC,1);
        out = hybZono(Gc,Gb,c,Ac,Ab,b);
end

end