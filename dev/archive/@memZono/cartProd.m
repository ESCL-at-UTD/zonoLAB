function obj = cartProd(obj1,obj2,cdsKeyPrefix)
    arguments
        obj1 memZono
        obj2 memZono
        cdsKeyPrefix = 'cds';
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
        obj1.c(idxd1);
        obj2.c(idxd2)
        ];
    A_ = [
        obj1.A(idxc1,idxk1), obj1.A(idxc1,idxks1), zeros(length(c1),length(k2));
        zeros(length(c2),length(k1)), obj2.A(idxc2,idxks2), obj2.A(idxc2,idxk2);
        ];
    b_ = [
        obj1.b(idxc1);
        obj2.b(idxc2);
        ];

    % hybrid Zono
    if obj1.vset(idxks1) ~= obj1.vset(idxks2)
        error('c/d factors not lining up');
        % TODO???: add option for c/d not matching => additonal
        % factors/constraints?
    end
    vset_ = [obj1.vset(idxk1),obj1.vset(idxks1),obj2.vset(idxk2)];

    % Labeling
    keys_.factors = [k1,ks,k2];
    keys_.dims = [d1,d2];
    keys_.cons = [c1,c2];

    %% Shared Dimensions
    if ~isempty(ds)
        % Intersection terms
        R = zeros(length(ds),obj1.n);
        for k = 1:length(ds)
            i = idxds2(k); j = idxds1(k); R(i,j) = 1; 
        end
        % cds key labels
        cds{length(ds)} = [];
        for i = 1:length(ds)
            cds{k} = sprintf('%s_%s_%d',cdsKeyPrefix,ds{k},k);
        end
        % Matrices
        G_ = [G_;
            obj1.G(idxds1,idxk1), obj1.G(idxds1,idxks1), zeros(length(ds),length(k2));
        ];
        c_ = [c_;
            obj1.c(idxds1);
        ];
        A_ = [A_; 
            R*obj1.G(:,idxk1), R*obj1.G(:,idxks1)-obj2.G(idxds2,idxks2), -obj2.G(idxds2,idxk2) 
        ];
        b_ = [b_;
            obj2.c(idxds2) - R*obj1.c;
        ];

        % Labels
        keys_.dims = [keys_.dims, ds];
        keys_.cons = [keys_.cons, cds];
    end

    %% Shared Constraints
    if ~isempty(cs)
        % @JONAS - needed to add the shared factors indicies.  Also needed to use 'all' since we're comparing arrays...
        if all(obj1.A(idxcs1,idxks1) ~= obj2.A(idxcs2,idxks2),'all') || all(obj1.b(idxcs1) ~= obj2.b(idxcs2),'all')
            error('Shared Constraints are not identical')
            % TODO???: add option for non-identical constraints
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
