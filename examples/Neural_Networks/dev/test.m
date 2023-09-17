
n0 = 1; % input dimension
d = 2; % input domain
X = hybZono(d*eye(n0),[],zeros([n0,1]),[],[],[]); % input hybZono set

L = [n0,3,2,1]; % layer sizes
Ws={}; % weight matrices
bs={}; % bias vectors
for i = 1:(length(L)-1)
    Ws{i} = rand([L(i+1),L(i)]);
    bs{i} = rand([L(i+1),1]);
end

[Z,XZ] = relu_neural_network(Ws,bs,100,X);

figure(1)
plot(XZ,'b',0.8);