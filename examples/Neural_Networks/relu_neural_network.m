function [Y,XY] = relu_neural_network(Ws,bs,a,X)

if nargin < 3
    a = 1000;
end

c = a/2*[0;1];
Gc = a/2*[1,1,0,0;1,0,0,0];
Gb = [0;0];
Ac = a/2*[1,1,-1,-1;1,0,-1,0];
Ab = a/2*[-2;-1];
b = a/2*[0;-1];

% Hybrid Zonotope of a single ReLU
relu = hybZono(Gc,Gb,c,Ac,Ab,b);

if nargin < 4
    [n2,n1] = size(Ws{1});
    X = hybZono(a*eye(n1),[],zeros([n1,1]),[],[],[]);
end

% Neural network hybrid zonotope
NN = X;

v_indices = {};
x_indices = {}; % x would index starting from zero, so this is shifted
x_indices{1} = 1:X.n;
idx_max = X.n;

for i = 1:(length(bs)-1)
    [n2,n1] = size(Ws{i});

    In2 = eye(n2);
    % [v1 x1 ... vn xn] ordering
    relu_layer_stacked = hybZono( ...
                          kron(In2,Gc), ...
                          kron(In2,Gb), ...
                          repmat(c,n2,1), ...
                          kron(In2,Ac), ...
                          kron(In2,Ab), ...
                          repmat(b,n2,1));

    T = [ kron(eye(n2),[1,0]) ; kron(eye(n2),[0,1]) ];
    % [v1 ... vn , x1 ... xn] ordering
    relu_layer = T*relu_layer_stacked;
    
    % Keep track of v and x indices for adding layer-to-layer connections below
    v_indices{i} = idx_max+(1:n2);
    idx_max = idx_max+n2;
    x_indices{i+1} = idx_max+(1:n2);
    idx_max = idx_max+n2;

    % cartesian product to add layer to the hybrid zonotope neural network
    NN = cartProd(NN,relu_layer);
end

% make layer-to-layer connections via constraints
for i = 1:(length(bs)-1)
    % v^i = W^{i-1} x^{i-1} + b*{i-1}
    Rx = selection_matrix(x_indices{i},NN); % recall the x index is shifted up one
    Rv = selection_matrix(v_indices{i},NN);
    NN = intersect_in_place(NN, Ws{i}*(Rx*NN)+bs{i} , Rv);
end

% make connections from final hidden layer to output
Rx = selection_matrix(x_indices{1},NN);
RxL = selection_matrix(x_indices{end},NN);

X = Rx*NN;  % this version of X has all the factors from the neural network
Y = Ws{end}*RxL*NN+bs{end};

XY = hybZono( ...
                [X.Gc;Y.Gc], ...
                [X.Gb;Y.Gb], ...
                [X.c;Y.c], ...
                X.Ac, ...
                X.Ab, ...
                X.b);

end

% This assumes the set of generators of X and Y are exactly the same!
function Z = intersect_in_place(X,Y,R)
	Z = X;
	
	if nargin < 3
		R = eye(size(Z.c,1));
	end
	
	Z.Ac = [ X.Ac  ; (R*X.Gc -Y.Gc)];
	Z.Ab = [ X.Ab ; (R*X.Gb -Y.Gb)];
	Z.b = [ X.b ; (Y.c - R*X.c) ];
end	

% Convenience function to build the selection matrices for intersection or projection
function [R] = selection_matrix(idx_range,Z)
    R = zeros([length(idx_range),Z.n]);
    R(:,idx_range) = eye(length(idx_range));
end
