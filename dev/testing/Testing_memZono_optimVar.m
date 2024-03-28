% Testing Optimization in memZono
clear

n = 2; ng = 2; nc = 1;
X_G = optimvar('X_g_',[n ng]);
X_c = optimvar('x_c_',[n 1]);
if nc>0
    X_A = optimvar('X_A_',[nc ng]);
    X_b = optimvar('X_b_',[nc 1]);
    X = conZono(X_G,X_c,X_A,X_b);
else
    X = zono(X_G,X_c);
end

X = memZono(X,'x');

%% Propogation
xLabel = @(k) sprintf('x_%d',k);
A = [1,2;3,4]; b = ones(n,1); N = 5;
X_{1} = memZono(X.Z,'x_1');
for k = 2:N
    X_{k} = X_{k-1}.linMap(A,xLabel(k)) + memZono(b,xLabel(k));
end
% store as a single one
X_all = vertcat(X_{:});

% select x_1
X_end = X_all(X_all.keysStartsWith(xLabel(N)).dimKeys)

%% Support Function

