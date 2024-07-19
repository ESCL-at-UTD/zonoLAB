%% Workspace Clear
clear; clc;

%% Fixed 2D Test Case
G1 = [1 0 1; 0 1 -1];
c1 = [0; 0];
z1 = zono(G1, c1);

disp('Zono2HPoly Computation Time (n = 2, nG = 3)');
tic;
[C1, d1] = zono2HPoly(z1);
toc;
HPoly_1 = Polyhedron('H', [C1 d1]);

figure; grid on;
plot(z1, 'blue', 0.2); hold on;
HPoly_1.plot('color', 'red'); alpha(0.2);
legend('Zonotope $z_1$','H-Poly $H_1$', 'interpreter', 'latex');

%% Random 2D Test Case
G2 = rand(2, 50);
c2 = rand(2, 1);
z2 = zono(G2, c2);

disp('Zono2HPoly Computation Time (n = 2, nG = 50)');
tic;
[C2, d2] = zono2HPoly(z2);
toc;
HPoly_2 = Polyhedron('H', [C2 d2]);

figure; grid on;
plot(z2, 'blue', 0.2); hold on;
HPoly_2.plot('color', 'red'); alpha(0.2);
legend('Zonotope $z_2$','H-Poly $H_2$', 'interpreter', 'latex', 'location', 'best');

%% Random 3D Test Case
nG = 20;
G3 = rand(3, nG);
c3 = [0; 0; 0];
z3 = zono(G3, c3);

str = sprintf('Zono2HPoly Computation Time (n = 3, nG = %d)', nG);
disp(str);
tic;
[C3, d3] = zono2HPoly(z3);
toc;
HPoly_3 = Polyhedron('H', [C3 d3]);
HPoly_3.minHRep;

disp('Plotting in 3D...')
figure; grid on;
plot(z3, 'blue', 0.2); hold on;
HPoly_3.plot('color', 'red'); alpha(0.2);
xlabel('$e_1$', 'interpreter', 'latex');
ylabel('$e_2$', 'interpreter', 'latex');
zlabel('$e_3$', 'interpreter', 'latex');
legend('Zonotope $z_3$','H-Poly $H_3$', 'interpreter', 'latex');
drawnow;

%Generate 2d projections
disp('Computing zonotope projections...')
 z3_12 = projection(z3, [1 2]);
 z3_13 = projection(z3, [1 3]);
 z3_23 = projection(z3, [2 3]);
 
disp('Computing polyhedron projection (1,2)...')
 HPoly_3_12 = HPoly_3.projection([1 2]);
 HPoly_3_12.minHRep;
 
disp('Computing polyhedron projection (1,3)...')
 HPoly_3_13 = HPoly_3.projection([1 3]);
 HPoly_3_13.minHRep;
 
disp('Computing polyhedron projection (2,3)...')
 HPoly_3_23 = HPoly_3.projection([2 3]);
 HPoly_3_23.minHRep;

disp('Plotting 2D projections...')
 figure; grid on;
 plot(z3_12, 'blue', 0.2); hold on;
 HPoly_3_12.plot('color', 'red'); alpha(0.2);
 legend('Zonotope $z_3$','H-Poly $H_3$', 'interpreter', 'latex');
 xlabel('$e_1$', 'interpreter', 'latex');
 ylabel('$e_2$', 'interpreter', 'latex');
 title('(1,2) Projection')
 drawnow;

figure; grid on;
plot(z3_13, 'blue', 0.2); hold on;
HPoly_3_13.plot('color', 'red'); alpha(0.2);
legend('Zonotope $z_3$','H-Poly $H_3$', 'interpreter', 'latex');
xlabel('$e_1$', 'interpreter', 'latex');
ylabel('$e_3$', 'interpreter', 'latex');
title('(1,3) Projection')
drawnow;

figure; grid on;
plot(z3_23, 'blue', 0.2); hold on;
HPoly_3_23.plot('color', 'red'); alpha(0.2);
legend('Zonotope $z_3$','H-Poly $H_3$', 'interpreter', 'latex');
xlabel('$e_2$', 'interpreter', 'latex');
ylabel('$e_3$', 'interpreter', 'latex');
title('(2,3) Projection')
drawnow;

%% Random 4D Test Case
clear;
nG = 200;
G3 = rand(4, nG);
c3 = zeros(4, 1);
z3 = zono(G3, c3);

str = sprintf('Zono2HPoly Computation Time (n = 4, nG = %d)', nG);
disp(str);
tic;
[C3, d3] = zono2HPoly(z3, true);
toc;

