function out = conZono2AHPoly(obj)

out = affineHPoly;
if isempty(obj.nC)
    out.c = obj.c;
    out.G = obj.G;
    out.A = [-eye(obj.nG);eye(obj.nG)];
    out.b = ones(2*obj.nG,1);
else
    T = null(obj.A,'r');  % Computes the right null space of A.
    s = pinv(obj.A)*obj.b;
    P = Polyhedron('H',[T ones(size(s,1),1)-s; -T ones(size(s,1),1)+s]); % Inequality constraints on \xi variables.
    P = minHRep(P); % Removes redundant halfspaces
    out.c = obj.c+obj.G*s;
    out.G = obj.G*T;
    out.A = P.H(:,1:end-1);
    out.b = P.H(:,end);
end
end


