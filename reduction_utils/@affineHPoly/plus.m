function out = plus(obj1,obj2)
if obj1.n ~= obj2.n
    disp(['Cannot add a AH-polytopes in ',num2str(obj1.n),...
        ' dimensions to a AH-polytopes in ',num2str(obj2.n),' dimensions!'])
else
    out = affineHPoly;
    out.c = obj1.c + obj2.c;
    out.G = [obj1.G obj2.G];
    out.A = blkdiag(obj1.A,obj2.A);
    out.b = [obj1.b; obj2.b];
end
end