% Testing Optimization in memZono
clear

n = 2; ng = 2; nc = 0;
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
X_end = X_all(X_all.keysStartsWith(xLabel(N)).dimKeys);


%% Support Function
X_all.supportFunc({'x_1_1'});

%% Transient Peak Calc
clear

n = 2; m = 1; p = n;
A = optimvar('A',[n n]);
B = optimvar('B',[n m]);
C = optimvar('C',[p n]);
% D = 0; %<=== ignore

% IC
X0_G = 2*eye(n);
X0_c = ones(n,1);
X0 = zono(X0_G,X0_c);
X_0 = memZono(X0,'x_0');

U0_G = 0.1*[eye(m), ones(m)];
U0_c = zeros(m,1);
U0 = zono(U0_G,U0_c);
U_0 = memZono(U0,'u_0');


% Propogation example/test
xLabel = @(k) sprintf('x_%d',k);
uLabel = @(k) sprintf('u_%d',k);
yLabel = @(k) sprintf('y_%d',k);

N = 5;
X_{1} = memZono(X0,xLabel(0)).linMap(A,xLabel(1)) ...
    + memZono(U0,uLabel(1)).linMap(B,xLabel(1));
for k = 1:N-1
    U_{k} = memZono(U0,uLabel(k));
    Y_{k} = X_{k}.linMap(C,yLabel(k));
    X_{k+1} = X_{k}.linMap(A,xLabel(k+1)) ...
        + U_{k}.linMap(B,xLabel(k+1));
end
X_all = vertcat(X_{:},Y_{:},U_{:});


% Support functions based on X_all
h_var_k = @(var,k,d) X_all.supportFunc( ...
    X_all.keysStartsWith(sprintf('%s_%d',var,k)).dimKeys, d);

C = eye(p,n); %<=== will be selecting frist row... so it's really grabbing just first state
C = C(1,:);
h_X_all_k = @(k,C) h_var_k(xLabel(k),C');


% Support Functions based on propogation independently
h_X_k = @(k,C) X_0.supportFunc('x_0',(C*A^k)');

h_X_k(1,C)




