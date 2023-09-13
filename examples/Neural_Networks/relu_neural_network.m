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
relu = hybZono(c,Gc,Gb,Ac,Ab,b);

if nargin < 4
    [n2,n1] = size(Ws{1});
    X = hybZono( zeros([n1,1]), a*eye(n1), [],[],[],[]);
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
    relu_layer_stacked = hybZono( repmat(c,n2,1), ...
                          kron(In2,Gc), ...
                          kron(In2,Gb), ...
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
    NN = stack_zonotopes(NN,relu_layer);
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
RvL = selection_matrix(v_indices{end},NN);

X = Rx*NN;  % this version of X has all the factors from the neural network
Y = Ws{end}*RvL*NN+bs{end};

XY = hybZono( [X.c;Y.c], ...
                [X.Gc;Y.Gc], ...
                [X.Gb;Y.Gb], ...
                X.Ac, ...
                X.Ab, ...
                X.b);

end

% This assumes the set of generators of X and Y are exactly the same!
function Z = intersect_in_place(X,Y,R)
    X = X.getDimensions;
    Y = Y.getDimensions;
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

% same as Cartesian product
% for some reason the OG hybZono toolbox wasn't doing the same thing?
function [Z] = stack_zonotopes(X,Y)
% Stacking two independent zonotopes. 
nx = max( size(X.Gc,1) , size(X.Gb,1) );
ny = max( size(Y.Gc,1) , size(Y.Gb,1) );
ncx = max( size(X.Ac,1) , size(X.Ab,1) );
ncy = max( size(Y.Ac,1) , size(Y.Ab,1) );

Z = hybZono( [X.c ; Y.c], ...
             [X.Gc , zeros(nx,size(Y.Gc,2)) ; zeros(ny,size(X.Gc,2)) , Y.Gc, ], ...
             [X.Gb , zeros(nx,size(Y.Gb,2)) ; zeros(ny,size(X.Gb,2)) , Y.Gb, ], ...
             [X.Ac , zeros(ncx,size(Y.Gc,2)) ; zeros(ncy,size(X.Gc,2)) , Y.Ac], ...
             [X.Ab , zeros(ncx,size(Y.Gb,2)) ; zeros(ncy,size(X.Gb,2)) , Y.Ab], ...
             [X.b ; Y.b]);
end
