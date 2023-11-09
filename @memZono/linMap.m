function obj = linMap(M,obj)
    arguments
        M double
        obj memZono
    end

    obj.G = M*obj.G;
    obj.c = M*obj.c;
end