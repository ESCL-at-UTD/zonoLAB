function out = plus(obj1,obj2)

% Standardized header

if ~strcmp(class(obj1),class(obj2)) % If classes do not match
    types = [isa(obj1,'hybZono') isa(obj2,'hybZono');...
             isa(obj1,'conZono') isa(obj2,'conZono');...
             isa(obj1,'zono')    isa(obj2,'zono')];
    if sum(types(1,:)) ~=0 % One of the sets is a hybrid zonotope
        if types(1,1) == 1
            out = plus(obj1,hybZono(obj2));
        else
            out = plus(hybZono(obj1),obj2);
        end
    elseif sum(types(2,:)) ~=0 % One of the sets is a constrained zonotope
        if types(2,1) == 1
            out = plus(obj1,conZono(obj2));
        else
            out = plus(conZono(obj1),obj2);
        end
    elseif sum(types(3,:)) ~=0 % One of the sets is a zonotope
        if types(3,1) == 1
            out = plus(obj1,zono(obj2));
        else
            out = plus(zono(obj1),obj2);
        end
    else
        warning(['Only hybrid zonotopes, constrained zonotopes, and zonotopes';...
                  'can be added with each other or shifted by vectors.'])
    end
else % Once classes match then sum
    switch class(obj1)
        case 'zono'
            G = [obj1.G obj2.G];
            c = obj1.c + obj2.c;
            out = zono(G,c);
        case 'conZono'
            G = [obj1.G obj2.G];
            c = obj1.c + obj2.c;
            A = blkdiag(obj1.A,obj2.A);
            b = [obj1.b; obj2.b];
            out = conZono(G,c,A,b);
        case 'hybZono'
            Gc = [obj1.Gc obj2.Gc];
            Gb = [obj1.Gb obj2.Gb];
            c = obj1.c + obj2.c;
            Ac = blkdiag(obj1.Ac,obj2.Ac);
            Ab = blkdiag(obj1.Ab,obj2.Ab);
            b = [obj1.b; obj2.b];
            out = hybZono(Gc,Gb,c,Ac,Ab,b);
    end
end
end