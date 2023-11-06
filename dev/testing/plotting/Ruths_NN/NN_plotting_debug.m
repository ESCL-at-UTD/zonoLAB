load('NN_plotting_debug.mat');

%% Plot training data
az = -138; el = 37;
figure(3)
subplot(3,2,1)
cla;
surf(TH,V,dx,'EdgeColor','none');
view(az,el);
xlabel('\theta')
ylabel('v')
zlabel('x')

subplot(3,2,2)
cla;
surf(TH,V,dy,'EdgeColor','none');
view(az,el);
xlabel('\theta')
ylabel('v')
zlabel('y')

%% Plot test data (the predictions of the neural network)
figure(3)
subplot(3,2,3)
cla;
surf(TH,V,dx_test,'EdgeColor','none');
view(az,el);
xlabel('\theta')
ylabel('v')
zlabel('x')

subplot(3,2,4)
cla;
surf(TH,V,dy_test,'EdgeColor','none');
view(az,el);
xlabel('\theta')
ylabel('v')
zlabel('y')

%% Plot the HybZono   <---- This is what is not working
figure(3)
subplot(3,2,5);
cla;
XZ_proj = Rxz_thvx*XZ;
R = rref([XZ_proj.Ac XZ_proj.Ab XZ_proj.b]);
XZ_proj.Ac = R(:,1:XZ_proj.nGc);
XZ_proj.Ab = R(:,1+XZ_proj.nGc:end-1);
XZ_proj.b = R(:,end);
plot(XZ_proj,'b',0.7);
hold on;
hold off;
grid on;
view(az,el);
xlabel('\theta'); ylabel('v'); zlabel('dx')
%%
subplot(3,2,6);
cla;
plot(Rxz_thvy*XZ,'b',0.7);
hold on;
hold off;
grid on;
view(az,el);
xlabel('\theta'); ylabel('v'); zlabel('dy')

%% Plot the HybZono using lines
figure(3)
subplot(3,2,5);
cla;
Rxz_vx = zeros(2,6); Rxz_vx(:,[2,4]) = eye(2);
for th=theta
theta_zono = hybZono([],[],[th],[],[],[]);
theta_strip = Rxz_vx*and(XZ,theta_zono,[1,0,0,0,0,0]);
theta_strip.Gc = [zeros(1,theta_strip.nGc);theta_strip.Gc];
theta_strip.Gb = [zeros(1,theta_strip.nGb);theta_strip.Gb];
theta_strip.c = [th;theta_strip.c];
plot(theta_strip,'r',1);
hold on;
end
hold off;
grid on;
view(az,el);
xlabel('\theta'); ylabel('v'); zlabel('dx')

subplot(3,2,6);
Rxz_vy = zeros(2,6); Rxz_vy(:,[2,5]) = eye(2);
for th=theta
theta_zono = hybZono([],[],[th],[],[],[]);
theta_strip = Rxz_vy*and(XZ,theta_zono,[1,0,0,0,0,0]);
theta_strip.Gc = [zeros(1,theta_strip.nGc);theta_strip.Gc];
theta_strip.Gb = [zeros(1,theta_strip.nGb);theta_strip.Gb];
theta_strip.c = [th;theta_strip.c];
plot(theta_strip,'r',1);
hold on;
end
hold off;
grid on;
view(az,el);
xlabel('\theta'); ylabel('v'); zlabel('dy')