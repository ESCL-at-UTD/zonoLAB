%% Testing Optimization integration
clear

n = 2; ng = 2;
X_G = optimvar('X_g_',n,ng);
X_c = optimvar('x_c_',n,1);

X = conZono(X_G,X_c)

A = [1,2;3,4];

X_{1} = X;
for k = 2:5
    X_{k} = A*X_{k-1};
end

X_{end}.G, X_{end}.c

H = eye(n); f = ones(n,1);
% X_half = halfspaceIntersection(X,H,f) %<== fails for optimvar() due to abs...

