V = [ -1.7 -0.4; -0.4  0.7; 1.2 -0.8; 0 0.8; 1.3 0.9; -0.3 0.6];
P4 = Polyhedron(V);
figure; hold on
plot(P4)

P4.H

H = P4.H(:,1:end-1);
f = P4.H(:,end);

Zc = hPoly2conZono([H f]);

plot(Zc,'b',0.5)

%%

% cZ = conZonotope(Zc.c,Zc.G,Zc.A,Zc.b)
% 
% figure;plot(cZ)