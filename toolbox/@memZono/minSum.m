function obj = minSum(obj1,obj2)
    arguments
        obj1 memZono
        obj2 memZono
    end
    
    % shared keys
    [k1,ks,k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);

    % Boolean maps for shared factors <--- make clearer?
    [~,idxk1] = ismember(k1,obj1.factorKeys);
    [~,idxks1] = ismember(ks,obj1.factorKeys);
    [~,idxks2] = ismember(ks,obj2.factorKeys);
    [~,idxk2] = ismember(k2,obj2.factorKeys);

    % Data selection
    G_ = [obj1.G(:,idxk1),obj1.G(:,idxks1)+obj2.G(:,idxks2), obj2.G(:,idxk2)];
    c_ = obj1.c + obj2.c;
    A_ = [obj1.A(:,idxk1), obj1.A(:,idxks1), zeros(obj1.nC,obj2.nG);
            zeros(obj2.nC,obj1.nG), obj2.A(:,idxks2), obj2.A(:,idxk2)];
    b_ = [obj1.b; obj2.b];
    vset_ = [obj1.vset(idxk1),obj1.vset(idxks1),obj2.vset(idxk2)]; %<--- add check for c/d the same?

    % Labeling
    keys_.factors = [k1,ks,k2];
    keys_.dims = obj1.dimKeys; %<--- add check for same?
    keys_.cons = [obj1.conKeys; obj2.conKeys];

    % Construct object
    obj = memZono(G_,c_,A_,b_,vset_,keys_);

end
