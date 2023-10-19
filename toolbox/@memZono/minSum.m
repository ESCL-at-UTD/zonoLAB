function obj = minSum(obj1,obj2)
    arguments
        obj1 memZono
        obj2 memZono
    end
    
    % shared keys
    [k1,ks,k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);

    % Boolean maps for shared factors <--- make clearer?
    bs1 = ismember(obj1.factorKeys,ks);
    bs2 = ismember(obj2.factorKeys,ks);

    % Data selection
    G_ = [obj1.G(:,~bs1),obj1.G(:,bs1)+obj2.G(:,bs2), obj2.G(:,~bs2)];
    c_ = obj1.c + obj2.c;
    A_ = [obj1.A(:,~bs1), obj1.A(:,bs1), zeros(obj1.nC,obj2.nG);
            zeros(obj2.nC,obj1.nG), obj2.A(:,bs2), obj2.A(:,~bs2)];
    b_ = [obj1.b; obj2.b];
    vset_ = [obj1.vset(~bs1),obj1.vset(bs1),obj2.vset(~bs2)]; %<--- add check for c/d the same?

    % Labeling
    keys_.factors = [k1,ks,k2];
    keys_.dims = obj1.dimKeys; %<--- add check for same?
    keys_.cons = [obj1.conKeys; obj2.conKeys];

    % Construct object
    obj = memZono(G_,c_,A_,b_,vset_,keys_);

end
