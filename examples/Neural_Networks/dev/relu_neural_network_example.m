clear all; close all; clc;

example = 'sincos';
%example = 'linear';

%% Create the Training/Testing Data

N = 200;

%Training Input
x1_train = linspace(-5,5,N);
x2_train = linspace(-5,5,N);
[X1_train,X2_train] = meshgrid(x1_train,x2_train); % input domain
[n,m] = size(X1_train);
input_train = [reshape(X1_train,n*m,1),reshape(X2_train,n*m,1)]; % reshape input data

% Testing input (more densely sampled)
x1_test = linspace(min(x1_train),max(x1_train),2*N);
x2_test = x1_test;
[X1_test,X2_test] = meshgrid(x1_test,x2_test); % input domain
[p,q] = size(X1_test);
input_test = [reshape(X1_test,p*q,1),reshape(X2_test,p*q,1)];

% Function to learn
if example == 'sincos'
    f = @(X1,X2)cos(X1)+sin(X2);
else
    f = @(X1,X2)3*X1;
end

% Training data and validation data
Y_train = f(X1_train,X2_train);
Y_validate = f(X1_test,X2_test);

%%
figure('Position',[0,0,1500,400])
subplot(1,4,1)
surf(X1_train, X2_train, Y_train, 'EdgeColor', 'none');
title('Actual Function')

%% Train Network
output_train = reshape(Y_train,n*m,1); % reshape output data  

if example == 'sincos'
    h_layers = [20,10,10];
else
    h_layers = [4,5,4];
end

layers = [featureInputLayer(2) 
         fullyConnectedLayer(h_layers(1))
         reluLayer
         fullyConnectedLayer(h_layers(2))
         reluLayer
         fullyConnectedLayer(h_layers(3))
         reluLayer
         fullyConnectedLayer(1)
         regressionLayer];
     
options = trainingOptions("sgdm",'Momentum', 0.95, MaxEpochs = 200,Plots='training-progress');
%net = trainNetwork(input_train, output_train, layers, options);

% .... or load it from a mat file
if example == 'sincos'
    load("sin_cos_2_20_10_10_1.mat");
else
    load("linear_3x1_2_4_5_4_1.mat");
end

%% Validate the Network
output_test = predict(net,input_test); 
Y_test = reshape(output_test,p,q);

subplot(1,4,2)
%plot3(X_test(:,1),X_test(:,2),output_test,'.')
surf(X1_test,X2_test,Y_test, 'EdgeColor', 'none');
title('Neural Net Approx')

subplot(1,4,3)
surf(X1_test,X2_test,Y_test-Y_validate, 'EdgeColor', 'none');
title('Approximation Error')

pause(0.1)

%% Plot neural network output space 

[x1_min, x1_max] = deal(double(min(x1_test)), double(max(x1_test)));
[x2_min, x2_max] = deal(double(min(x2_test)), double(max(x2_test)));
domain = [x1_min,x1_max,x2_min,x2_max];
x1 = (x1_max - x1_min)*rand(1, N) + x1_min;
x2 = (x2_max - x2_min)*rand(1, N) + x2_min;

Z_nn = [x1; x2];

Ws = {double(net.Layers(2).Weights), double(net.Layers(4).Weights), double(net.Layers(6).Weights), double(net.Layers(8).Weights)};
bs = {double(net.Layers(2).Bias), double(net.Layers(4).Bias), double(net.Layers(6).Bias), double(net.Layers(8).Bias)};

%% Construct Hybrid Zonotope
tic

% Construct input zonotope
g11 = (x1_max - x1_min)/2;
g22 = (x2_max - x2_min)/2;
Gx = diag([g11, g22]);
cx = zeros(2, 1);
X = hybZono(Gx, [], cx, [], [], []);

if example == 'sincos'
    a = 1000;
else
    a = 100;
end
[Z,XZ] = relu_neural_network(Ws,bs,a,X);

fprintf('Zonotope model: ')
toc

%% Plot Hybrid Zonotope
subplot(1,4,4)
plot(XZ,'r',0.8);
grid on;
title('Hybrid Zonotope')
toc

%%
% if example == 'sincos'
%     save("sin_cos_2_20_10_10_1.mat");
%     save("sin_cos_2_20_10_10_1_W_b.mat","bs","Ws","domain","a");
% else
%     save("linear_3x1_2_4_5_4_1.mat");
%     save("linear_3x1_2_4_5_4_1_W_b.mat","bs","Ws","domain","a");
% end