function obj = cartProd(obj1,obj2)
    arguments
        obj1 memZono
        obj2 memZono
    end

    %% Keys
    % shared factors
    [k1,ks,k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);
    [idxk1,idxks1,idxks2,idxk2] = memZono.getKeyIndices(obj1.factorKeys,obj2.factorKeys);
    % shared dims
    if isempty(obj1.dimKeys); obj1.dimKeys = 'd1'; end
    if isempty(obj2.dimKeys); obj2.dimKeys = 'd2'; end
    [d1,ds,d2] = memZono.getUniqueKeys(obj1.dimKeys,obj2.dimKeys);
    [idxd1,idxds1,idxds2,idxd2] = memZono.getKeyIndices(obj1.dimKeys,obj2.dimKeys);
    % shared cons
    [c1,cs,c2] = memZono.getUniqueKeys(obj1.conKeys,obj2.conKeys);
    [idxc1,idxcs1,idxcs2,idxc2] = memZono.getKeyIndices(obj1.conKeys,obj2.conKeys);
    % R for repeated dims
    R = zeros(length(ds),obj1.n);
    for k = 1:length(ds)
        i = idxds2(k); j = idxds1(k); R(i,j) = 1;
    end

    %% Data selection
    G_ = [
        obj1.G(idxd1,idxk1), obj1.G(idxd1,idxks1), zeros(length(d1),length(k2));
        obj1.G(idxds1,idxk1), obj1.G(idxds1,idxks1), zeros(length(ds),length(k2));
        zeros(length(d2),length(k1)), obj2.G(idxd2,idxks2), obj2.G(idxd2,idxk2)
        ];
    c_ = [
        obj1.c(idxd1);
        obj1.c(idxds1);
        obj2.c(idxd2)
        ];
    A_ = [
        obj1.A(idxc1,idxk1), obj1.A(idxc1,idxks1), zeros(length(c1),length(k2));
        obj1.A(idxcs1,idxk1), obj1.A(idxcs1,idxks1), zeros(length(cs),length(k2));
        zeros(length(cs),length(k1)), obj2.A(idxcs2,idxks2), obj2.A(idxcs2,idxk2);
        zeros(length(c2),length(k1)), obj2.A(idxc2,idxks2), obj2.A(idxc2,idxk2);
        % R*obj1.G(:,idxk1), R*obj1.G(:,idxks1)-obj2.G(idxds2,idxks2), -obj2.G(idxds2,idxk2) %<-- intersection term
        ];
    b_ = [
        obj1.b(idxc1);
        obj1.b(idxcs1);
        obj2.b(idxcs2);
        obj2.b(idxc2);
        % obj2.c(idxds2) - R*obj1.c;
        ];
    if ~isempty(idxds2)
        A_ = [
            A_;
            R*obj1.G(:,idxk1), R*obj1.G(:,idxks1)-obj2.G(idxds2,idxks2), -obj2.G(idxds2,idxk2) %<-- intersection term
        ];
        b_ = [b_;
            obj2.c(idxds2) - R*obj1.c;
        ];
    end
    vset_ = [obj1.vset(idxk1),obj1.vset(idxks1),obj2.vset(idxk2)]; %<--- add check for c/d the same?

    %% Labeling
    keys_.factors = [k1,ks,k2];
    keys_.dims = [d1,ds,d2]; %<--- add check for same?
    % Constraints
    if isempty(cs)
        cs1 = cs; cs2 = [];
    else
        cs1 = cs; cs2{length(cs)} = [];
        for i = 1:length(cs)
            cs2{i} = sprintfs('%s_%d',cs{i},2);
        end
    end
    if isempty(ds)
        cds = [];
    else
        cds{length(ds)} = [];
        for i = 1:length(ds)
            max_iter=1000;
            for j = 1:max_iter
                cds{i} = sprintf('i%i_%s',j,ds{i});
                if all(~ismember(cds{i},c1)) && all(~ismember(cds{i},cs)) && all(~ismember(cds{i},c2))
                    cds{i} = sprintf('i%i_%s',j,ds{i});
                    break
                end
                if j == max_iter
                    fprintf('ERROR: constraint at max number of iterations\n');
                end
            end
        end
    end
    keys_.cons = [c1,cs1,cs2,c2,cds];

    % Define memZono
    obj = memZono(G_,c_,A_,b_,vset_,keys_);
    
end
