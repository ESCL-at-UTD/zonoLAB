clear all, close all, clc
rng_seed = 42;
rng(rng_seed);

X.c = zeros(2,1);
X.Gc = randn(2,4);
X.Gb = 10*randn(2,6);
X.Ac = [randn(1,4); zeros(1,4)];
X.Ab = [zeros(1,6); 1 0 1 0 0 0];
X.b = zeros(2,1);
X = hybZono(X.Gc, X.Gb, X.c, X.Ac, X.Ab, X.b);

Y.c = 10*randn(2,1)+[30;20];
Y.Gc = diag([1 1])*randn(2,4);
Y.Gb = 20*randn(2,2);
Y.Ac = rand(1,4);
Y.Ab = zeros(1,2);
Y.b = 0;
Y = hybZono(Y.Gc, Y.Gb, Y.c, Y.Ac, Y.Ab, Y.b);

Z = convexHull(X,Y);

figure
hold on, grid on, grid minor
plot(X, 'r', 0.5)
plot(Y, 'b', 0.5)
plot(Z, 'g', 0.3)
legend({'$\mathcal{X}$','$\mathcal{Y}$','$\mathcal{Z} = \textbf{hull}(\mathcal{X} \cup \mathcal{Y})$'}, 'Location','eastoutside','Interpreter','latex')

%%
clear, close all, clc

X = randomSet(1, 'conZono', 3, 4, 4, 1);
Y = randomSet(2, 'conZono', 3, 4, 4, 1);
Z = convexHull(X, Y);
figure, hold on, grid on, grid minor
plot(X, 'b', 0.5)
plot(Y, 'r', 0.5)
plot(Z, 'g', 0.3)

%%
clear, close all, clc
rng_seed = 44;
rng(rng_seed);

X.c = zeros(2,1);
X.Gc = randn(2,4)+[0;1];
X.Gb = 10*randn(2,5);
X.Ac = [randn(1,4); zeros(1,4)];
X.Ab = [zeros(1,5); 1 1 0 1 -1];
X.b = zeros(2,1);
X = hybZono(X.Gc, X.Gb, X.c, X.Ac, X.Ab, X.b);
Y = randomSet(44, 'conZono', 2, randi(10), [], randi(2));
Y.G = Y.G*10;
Y.c = Y.c + [20; -24];
Z = convexHull(X, Y);

figure
hold on, grid on, grid minor
plot(X, 'r', 0.7)
plot(Y, 'b', 0.7)
plot(Z, 'g', 0.2)
