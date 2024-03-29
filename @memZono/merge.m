% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       A dimension-aware and memory-enabled intersection of memZono objects
%   Syntax:
%       [Z] = merge(X,Y,sharedDimLabels)
%   Inputs:
%       X               - memZono in R^n
%       Y               - memZono in R^m
%       sharedDimLabels - either (1) a cell array of new constraint key 
%                                    labels
%                                (2) a string prefix to use for new 
%                                    constraint key labels
%   Outputs:
%       Z - memZono in R^p, where n,m < p <= n+m
%           Shared dimensions are intersected, unshared dimensions kept
%   Notes:
%       Shared dimensions undergo an intersection and unshared dimensions 
%       are maintained. Intersections will result in additional constraints 
%       provided and labeled acording to sharedDimLabels.
%       Factors are aligned to preserve memory.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function obj = merge(obj1,obj2,sharedDimLabels)
    arguments
        obj1
        obj2
        sharedDimLabels = {};
    end

    % Input Conditioning
    if ~isa(obj1,'memZono') || ~isa(obj2,'memZono')
        error('Inputs must both be memZono objects')
    end
    
    %% Keys
    % shared factors
    [k1,ks,k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);
    [idxk1,idxks1,idxks2,idxk2] = memZono.getKeyIndices(obj1.factorKeys,obj2.factorKeys);
    % shared dims
    [d1,ds,d2] = memZono.getUniqueKeys(obj1.dimKeys,obj2.dimKeys);
    [idxd1,idxds1,idxds2,idxd2] = memZono.getKeyIndices(obj1.dimKeys,obj2.dimKeys);
    % shared cons
    [c1,cs,c2] = memZono.getUniqueKeys(obj1.conKeys,obj2.conKeys);
    [idxc1,idxcs1,idxcs2,idxc2] = memZono.getKeyIndices(obj1.conKeys,obj2.conKeys);

    %% Factor-based Memory Cartisian Product
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

    % hybrid Zono
    if obj1.vset(idxks1) ~= obj2.vset(idxks2)
        error('c/d factors not lining up');
    end
    vset_ = [obj1.vset(idxk1),obj1.vset(idxks1),obj2.vset(idxk2)];

    % Labeling
    keys_.factors = [k1,ks,k2];
    keys_.dims = [d1,d2];
    keys_.cons = [c1,c2];

    %% Shared Dimensions
    if ~isempty(ds)
        % Intersection matrix definition
        R = zeros(length(ds),obj1.n);
        for k = 1:length(ds)
            j = idxds1(k); 
            R(k,j) = 1; 
        end
        
        % Interesecting Matrices
        G_ = [G_;
            obj1.G(idxds1,idxk1), obj1.G(idxds1,idxks1), zeros(length(ds),length(k2));
        ];
        c_ = [c_;
            obj1.c(idxds1,:);
        ];
        A_ = [A_; 
            R*obj1.G(:,idxk1), R*obj1.G(:,idxks1)-obj2.G(idxds2,idxks2), -obj2.G(idxds2,idxk2) 
        ];
        b_ = [b_;
            obj2.c(idxds2,:) - R*obj1.c;
        ];

        % Labels
        keys_.dims = [keys_.dims, ds];

        % Constraint Labels
        if ~isa(sharedDimLabels,'cell')
            if ~(isstring(sharedDimLabels)||ischar(sharedDimLabels))
                error('Intersection operation will add additional constraints but labels for new constraints are not given.'); 
            end
            cds{length(ds)} = [];
            for k = 1:length(ds)
                % TODO: Warning if these new conKey labels already exist in either obj1 or obj2
                cds{k} = sprintf('%s_%s_%d',sharedDimLabels,ds{k},k);
            end
        elseif length(sharedDimLabels) ~= length(ds)
            error('Intersection operation will add additional constraints but labels for new constraints are not given (or not enough given).');
        else
            cds = sharedDimLabels;
        end
        if ~isempty([obj1.conKeys,obj2.conKeys])
            if any(ismember(cds,[obj1.conKeys,obj2.conKeys]))
                error('Intersection labels already exists... provide new names')
            end
        end
        keys_.cons = [keys_.cons, cds];
    end

    %% Shared Constraints
    if ~isempty(cs)
        if all([isnumeric(obj1.A),isnumeric(obj2.A),isnumeric(obj1.b),isnumeric(obj2.b)])
            if all(obj1.A(idxcs1,idxks1) ~= obj2.A(idxcs2,idxks2),'all') || all(obj1.b(idxcs1) ~= obj2.b(idxcs2),'all')
                    error('Shared Constraints are not identical')
            end
        end
        A_ = [A_;
            zeros(length(cs),length(k1)), obj1.A(idxcs1,idxks1), zeros(length(cs),length(k2))
        ];
        b_ = [b_;
            obj1.b(idxcs1,:)
        ];
        % Labeling
        keys_.cons = [keys_.cons,cs];
    end

    %% Define memZono
    obj = memZono(G_,c_,A_,b_,vset_,keys_);
    
end
