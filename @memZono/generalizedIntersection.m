function obj = generalizedIntersection(obj1,obj2,R)
    arguments
        obj1 memZono
        obj2 memZono
        R = [];
    end

    if isempty(R); R = eye(obj2.n,obj1.n); end

    %% Keys
    % shared factors
    [k1,ks,k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);
    [idxk1,idxks1,idxks2,idxk2] = memZono.getKeyIndices(obj1.factorKeys,obj2.factorKeys);
    
    
    
    % TODO?????? add shared dims/cons?
    % % shared dims
    % if isempty(obj1.dimKeys); obj1.dimKeys = 'd1'; end
    % if isempty(obj2.dimKeys); obj2.dimKeys = 'd2'; end
    % [d1,ds,d2] = memZono.getUniqueKeys(obj1.dimKeys,obj2.dimKeys);
    % [idxd1,idxds1,idxds2,idxd2] = memZono.getKeyIndices(obj1.dimKeys,obj2.dimKeys);
    % % shared cons
    % [c1,cs,c2] = memZono.getUniqueKeys(obj1.conKeys,obj2.conKeys);
    % [idxc1,idxcs1,idxcs2,idxc2] = memZono.getKeyIndices(obj1.conKeys,obj2.conKeys);
    % R for repeated dims
    % R = zeros(length(ds),obj1.n);
    % for k = 1:length(ds)
    %     i = idxds2(k); j = idxds1(k); R(i,j) = 1;
    % end

    %% Data selection
    G_ = [obj1.G(:,idxk1), obj1.G(:,idxks1), zeros(obj1.n,length(k2))];
    c_ = obj1.c;
    A_ = [
        obj1.A(:,idxk1), obj1.A(:,idxks1), zeros(obj1.nC,length(k2));
        zeros(obj2.nC,length(idxk1)), obj2.A(:,idxks2), obj2.A(:,idxk2);
        R*obj1.G(:,idxk1), R*obj1.G(:,idxks1) - obj2.G(:,idxks2), -obj2.G(:,idxk2)
        ];
    b_ = [
        obj1.b;
        obj2.b;
        obj2.c - R*obj1.c;
        ];
    vset_ = [obj1.vset(idxk1),obj1.vset(idxks1),obj2.vset(idxk2)]; %<--- add check for c/d the same?

    % Labeling
    keys_.factors = [k1,ks,k2];
    keys_.dims = obj1.dimKeys; %<--- add check for same?
    % Constraints
    cons_new{size(R,1)} = [];
    for i = 1:size(R,1)
        cons_new{i} = sprintf('intersect_%d',i);
    end
    keys_.cons = [obj1.conKeys,obj2.conKeys,cons_new];

    % Define memZono
    obj = memZono(G_,c_,A_,b_,vset_,keys_);
    
end