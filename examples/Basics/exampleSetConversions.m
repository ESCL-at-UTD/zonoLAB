% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example:
%   The conversion of a H-rep polytope to a constrained zonotope, a
%   collection of H-rep polytopes to a hybrid zonotope, and a collection of
%   vertices to a hybrid zonotope
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%% H-rep to conZono
% Example set from: https://scaron.info/blog/polyhedra-and-polytopes.html 
H = [0  1;...
     5 -2;...
    -1 -3;...
    -4 -2];
f = [7 36 -14 -26]';

Zc = conZono([H f]);

figure;
plot(Zc,'b',0.1)

%% H-rep collection to hybZono
rng(1)
a = 0;
b = 2*pi;
nSets = 5;
angles = (b-a).*rand(nSets,1) + a;

figure;hold on
for i = 1:nSets
    rotatedH = H*[cos(angles(i)) -sin(angles(i)); sin(angles(i)) cos(angles(i))];
    Zc = conZono([rotatedH f]);
    plot(Zc,'r',0.1)
    H_collection{i,1} = [rotatedH f];
end

Zh = hPoly2hybZono(H_collection);

plot(Zh,'b',0.1)

%% Collection of vertices to hybZono

V = [0 1 0.9 0;0 0 0.9 1];
nV = size(V,2);

figure;hold on
for i = 1:size(V,2)
    plot(V(1,i),V(2,i),'ok')
end

% One convex set
M = ones(4,1);
Zh = vPoly2hybZono({V,M});
plot(Zh,'b',0.1)
Zc = conZono([Zh.Gc Zh.Gb],Zh.c,[Zh.Ac Zh.Ab],Zh.b);
plot(Zc,'g',0.5)

% % Boundary of the set
% M = zeros(nV);
% M(1,1) = 1; M(1,nV) = 1;
% for i = 2:nV
%     M(i,i-1) = 1; 
%     M(i,i) = 1;
% end
% Zh = vPoly2hybZono({V,M});
% plot(Zh,'r',1)

% % Individual vertices
% M = eye(4);
% Zh = vPoly2hybZono({V,M});
% plot(Zh,'k',1)






