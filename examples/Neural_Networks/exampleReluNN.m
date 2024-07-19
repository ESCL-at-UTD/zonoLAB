% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example:
%   Exactly representing a neural network composed entirely of ReLU 
%   activation functions as a hybrid zonotope.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

% Parameters
% Load parameters from mat file
% This neural network approximates the function f = cos(X1)+sin(X2) over the domain [-5,5]x[-5,5]
load("relu_sin_cos_2_20_10_10_1.mat");

% Plotting function over the trained dataset
figure('Position',[0,0,1500,300])
subplot(1,4,1)
surf(X1_train, X2_train, Y_train, 'EdgeColor', 'none');
title('Actual Function')
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
zlabel('$y$','interpreter','latex')
zlim([-2 2])
set(gca,'fontsize',18,'fontname','times new roman')

% Plotting function over the test dataset
subplot(1,4,2)
surf(X1_test,X2_test,Y_test, 'EdgeColor', 'none');
title('Neural Net Approx')
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
zlabel('$\hat{y}$','interpreter','latex')
zlim([-2 2])
set(gca,'fontsize',18,'fontname','times new roman')

% Plotting the approximation error between the actual function and the trained neural network
subplot(1,4,3)
surf(X1_test,X2_test,Y_test-Y_validate, 'EdgeColor', 'none');
title('Approximation Error')
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
zlabel('$y-\hat{y}$','interpreter','latex')
zlim([-2 2])
set(gca,'fontsize',18,'fontname','times new roman')

% Construct input zonotope
cdomain = num2cell(domain);
[x1_min, x1_max, x2_min, x2_max] = cdomain{:};
g11 = (x1_max - x1_min)/2;
g22 = (x2_max - x2_min)/2;
Gx = diag([g11, g22]);
cx = zeros(2, 1);
X = hybZono(Gx, [], cx, [], [], []);

% Construct the output zonotope Z and the lifted input-output mapping XZ
fprintf('Zonotope model: ')
tic
[NN,Y] = reluNN(X,Ws,bs,a);
toc

% Plot Hybrid Zonotope
subplot(1,4,4)
% figure
plot(NN,'r',1);
grid on;
title('Hybrid Zonotope')
xlabel('$x_1$','interpreter','latex')
ylabel('$x_2$','interpreter','latex')
zlabel('$\hat{y}$','interpreter','latex')
zlim([-2 2])
set(gca,'fontsize',18,'fontname','times new roman')
toc

% exportgraphics(gcf,'cosSinNN.pdf','ContentType','vector')