function obj = linMap(obj,M)
    arguments
        obj memZono
        M
    end

    obj.G = M*obj.G;
    obj.c = M*obj.c;
end