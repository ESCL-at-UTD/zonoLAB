% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Construct a hybZono representing the same set as the input, but
%       whose convex relaxation is equal to the convex hull (the
%       tightest possible convex relaxation). This property is known in the
%       literature as the representation being "sharp".
%
%       This is an implementation of the method for constructing sharp
%       representations as given in:
%       Sherali, H. D., ; Adams, W. P. (1994). A hierarchy of relaxations 
%       and convex hull characterizations for mixed-integer zero-one 
%       programming problems. Discrete Applied Mathematics, 52, 83–106.
%       https://www.sciencedirect.com/science/article/pii/0166218X9200190W
%   Syntax:
%       Z = sharpHybZono(X)
%   Inputs:
%       X   - hybZono object
%   Outputs:
%       Z   - hybZono object that represents the same set as X, and has
%       the property that its convex relaxation is its convex hull
%   Notes:
%       The resulting hybrid zonotope will have the following complexity:
%           Z.nGb = X.nGb
%           Z.nGc = 2^X.nGb * (3*X.nGc + 1) - 1 - X.nGb
%           Z.nC  = 2^X.nGb * ( X.nC + 2*X.nGc)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function Z = sharpHybZono(X)

    % convert the original matrices to a {0,1}-HZ instead of {-1,1}-HZ
    Abar = 2*X.Ab;
    Bbar = 2*X.Ac;
    bbar = X.b + X.Ac*ones(X.nGc, 1) + X.Ab*ones(X.nGb, 1);
    
    % determine key dimensions and parameters from the data
    n = size(Abar, 2);
    m = size(Bbar, 2);
    d = n;
    D = min([n, d+1]);
    R = size(Abar, 1);

    N = 1:n;
    allJ = PowerSet(N);

    % total number of constraints that will be present in our hull set
    nC = R*length(allJ) + 2*m*(2^d*nchoosek(n,d));

    Aw = zeros(nC, 2^n-1);
    Av = zeros(nC, m*(2^n));
    As = zeros(nC, 2*m*(2^d*nchoosek(n, d)));
    b = zeros(nC, 1);

    % for all J that are subsets (order <= d) of N
    for i = 1:length(allJ)
        J = allJ{i};
        % coefficient on w_J is sum(alpha_rj)-b_r
        if ~isempty(J)
            Aw((i-1)*R+1:R*i, i-1) = sum(Abar(:, J), 2) - bbar;
        else
            b((i-1)*R+1:R*i) = bbar;
        end
        notJ = setdiff(N, J);
        for j = notJ
            % coefficient on w_Juj
            Aw((i-1)*R+1:R*i, getWDims([J j], allJ)) = Abar(:, j);
        end
        
        % coefficient on v_Jk
        Av((i-1)*R+1:R*i, ((i-1)*m+(1:m))) = Bbar;
    end

    [allJ1, allJ2] = makeJ1J2(n, d);
    idx = R*length(allJ)+1;
    s = 1;
    % for all (J1, J2) of order d
    for i = 1:length(allJ1)
        J1 = allJ1{i};
        J2 = allJ2{i};
        J2subsets = PowerSet(J2);
        
        for k = 1:m
            if ~isempty(J1)
                for ii = 1:length(J2subsets)
                    I = J2subsets{ii};
                    % f_d(J1, J2) - f_d^k(J1, J2) = slack variable
                    Aw(idx, getWDims([J1 I], allJ)) = (-1)^numel(I);
                    Av(idx, (getWDims([J1 I], allJ)*m+k)) = -(-1)^numel(I);

                    % f_d^k(J1, J2) = slack variable
                    Av(idx+1, (getWDims([J1 I], allJ)*m+k)) = (-1)^numel(I);
                end
                
            else
                b(idx) = -1;
                Av(idx, k) = -1;
                Av(idx+1, k) = 1;
                for ii = 2:length(J2subsets) % start at 2 to ignore empty set
                    I = J2subsets{ii};
                    Aw(idx, getWDims(I, allJ)) = (-1)^numel(I);
                    Av(idx, (getWDims(I, allJ)*m+k)) = -(-1)^numel(I);  
                    Av(idx+1, (getWDims([J1 I], allJ)*m+k)) = (-1)^numel(I);
                end
            end
            As(idx, s) = -1;
            As(idx+1, s+1) = -1;
            s = s+2;
            idx = idx+2;
        end
    end

    A_collect = [Aw Av As];
    numVars = size(A_collect, 2);
    Gw = [X.Gb zeros(X.n, size(Aw, 2)-size(X.Gb, 2))];
    Gv = [X.Gc zeros(X.n, size(Av, 2)-size(X.Gc, 2))];
    Gs = zeros(X.n, size(As, 2));
    Anew = 0.5*A_collect;
    bnew = b - .5*A_collect*ones(numVars, 1);

    % According to Thm3.2 in Sherali, if you keep all of the original
    % binaries as binary, then you get a (sharp) representation of the 
    % original set. If you relax them all to continuous, then you get the
    % convex hull. 
    Z = hybZono([Gw(:, n+1:end) Gv Gs], Gw(:, 1:n), X.c, Anew(:, n+1:end), Anew(:,1:n), bnew);    
end

function [ P ] = PowerSet( S )
    n = numel(S);
    x = 1:n;
    P = cell(1,2^n);
    p_ix = 2;
    for nn = 1:n
        a = combnk(x,nn);
        for j=1:size(a,1)
            P{p_ix} = S(a(j,:));
            p_ix = p_ix + 1;
        end
    end
end

function out = getWDims(idxs, allJ)
    idxs = sort(idxs);
    for i = 1:length(allJ)
        if isequal(allJ{i}, idxs)
            out = i-1;
            return
        end
    end
    out = [];
end

function [J1, J2] = makeJ1J2(n, d)
    OPTS = nchoosek(1:n, d);
    idx = 0;
    bins = dec2bin(0:2^d-1)-'0';
    for i = 1:nchoosek(n,d)
        for j = 1:2^d
            idx = idx+1;
            J1{idx} = [];
            J2{idx} = [];
            thisbin = bins(j,:);
            for k = 1:length(thisbin)
                if thisbin(k)==0
                    J1{idx}(end+1) = OPTS(i, k);
                else
                    J2{idx}(end+1) = OPTS(i, k);
                end
            end
        end
    end
end