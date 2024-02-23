% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       A dimension-aware and memory-preserving minkowski sum of memZono objects
%   Syntax:
%       [Z] = combine(X,Y)
%   Inputs:
%       X - memZono in R^n
%       Y - memZono in R^m
%   Outputs:
%       Z - memZono in R^p, where n,m < p <= n+m
%           shared dimensions are summed, unshared dimensions are kept
%   Notes:
%       Shared dimensions undergo a minkowski sum operation while 
%       all other dimensions are kept. 
%       Factors are aligned to preserve memory.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function obj = combine(obj1,obj2)

    %% Keys
    % shared factors
    [k1,ks,k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);
    [idxk1,idxks1,idxks2,idxk2] = memZono.getKeyIndices(obj1.factorKeys,obj2.factorKeys);
    % shared dimensions
    [d1,ds,d2] = memZono.getUniqueKeys(obj1.dimKeys,obj2.dimKeys);
    [idxd1,idxds1,idxds2,idxd2] = memZono.getKeyIndices(obj1.dimKeys,obj2.dimKeys);
    % shared constraints
    [c1,cs,c2] = memZono.getUniqueKeys(obj1.conKeys,obj2.conKeys);
    [idxc1,idxcs1,idxcs2,idxc2] = memZono.getKeyIndices(obj1.conKeys,obj2.conKeys);

    %% Factor-based Memory Cartesian Product
    G_ = [
        obj1.G(idxd1,idxk1), obj1.G(idxd1,idxks1), zeros(length(d1),length(k2));
        zeros(length(d2),length(k1)), obj2.G(idxd2,idxks2), obj2.G(idxd2,idxk2)
        ];
    c_ = [
        obj1.c(idxd1,:);
        obj2.c(idxd2,:)
        ];
    A_ = [
        obj1.A(idxc1,idxk1), obj1.A(idxc1,idxks1), zeros(length(c1),length(k2));
        zeros(length(c2),length(k1)), obj2.A(idxc2,idxks2), obj2.A(idxc2,idxk2);
        ];
    b_ = [
        obj1.b(idxc1,:);
        obj2.b(idxc2,:);
        ];

    if obj1.vset(idxks1) ~= obj2.vset(idxks2)
        error('c/d factors not lining up');
    end
    vset_ = [obj1.vset(idxk1),obj1.vset(idxks1),obj2.vset(idxk2)];

    % Labeling
    keys_.factors = [k1,ks,k2];
    keys_.dims = [d1,d2];
    keys_.cons = [c1,c2];

    %% Shared Dimensions 
    % Minkowski Sum (modified for shared factors) is applied to shared dimensions
    if ~isempty(ds)
        % Matrices
        G_ = [G_;
            obj1.G(idxds1,idxk1), obj1.G(idxds1,idxks1)+obj2.G(idxds2,idxks2), obj2.G(idxds2,idxk2);
        ];
        c_ = [c_;
            obj1.c(idxds1)+obj2.c(idxds2);
        ];
        % Constraints do not change with Minkowski Sum

        % Labeling
        keys_.dims = [keys_.dims, ds];
    end

    %% Shared Constraints
    if ~isempty(cs)
        % TODO: The following if seems fragile.  Not clear constraints will always be exactly the same coefficients? (could be multiplied by a scalar?)
        if all(obj1.A(idxcs1,idxks1) ~= obj2.A(idxcs2,idxks2),'all') || all(obj1.b(idxcs1) ~= obj2.b(idxcs2),'all')
            error('Shared Constraints are not identical')
        end
        A_ = [A_;
            zeros(length(cs),length(k1)), obj1.A(idxcs1,idxks1), zeros(length(cs),length(k2))
        ];
        b_ = [b_;
            obj1.b(idxcs1)
        ];

        % Labeling
        keys_.cons = [keys_.cons,cs];
    end

    %% Define memZono
    obj = memZono(G_,c_,A_,b_,vset_,keys_);
end
