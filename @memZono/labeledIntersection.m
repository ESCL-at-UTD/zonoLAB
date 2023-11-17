function obj = labeledIntersection(obj1,obj2,dims)
    arguments
        obj1 memZono
        obj2 memZono
        dims = [];
    end

    % Select Dims if not defined....
    if isempty(obj1.dimKeys); obj1.dimKeys = 'd1'; end
    if isempty(obj2.dimKeys); obj2.dimKeys = 'd2'; end
    if isempty(dims)
        [~,dims,~] = memZono.getUniqueKeys(obj1.dimKeys,obj2.dimKeys);
    end
    [~,idxd1] = ismember(dims,obj1.dimKeys);
    [~,idxd2] = ismember(dims,obj2.dimKeys);
    if isempty(idxd1) || isempty(idxd2)
        error('No matching dimensions and/or dims unspecified')
    end

    % Generate R
    R = zeros(length(dims),obj1.n);
    for k = 1:length(dims)
        i = idxd2(k); j = idxd1(k); R(i,j) = 1;
    end

    obj = generalizedIntersection(obj1,obj2.projection(dims),R);
end