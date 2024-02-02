% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the hybZono representation of a ReLU neural network (NN).
%   Syntax:
%       [NN]   = reluNN(X,Ws,bs,a)
%       [NN,Y] = reluNN(X,Ws,bs,a)
%   Inputs:
%       X  - zonotopic set in R^n (hybZono, conZono, or zono object)
%       Ws - 1 x nL cell array of weight matrices for the nL layers
%       bs - 1 x nL cell array of bias vectors for the nL layers
%       a  - positive scalar bounding the domain of each activation
%            function to the interval [-a,a]
%   Outputs:
%       NN - hybZono in R^m representing NN function
%       Y  - hybZono in R^m representing output set of NN function
%            corresponding to input set X
%   Notes:
%       Assumes all activation functions are ReLU. Assumes the input to
%       each activation function is within the interval [-a,a].
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function varargout = reluNN(X,Ws,bs,a)

nOutputs = nargout;
if (nOutputs == 1) || (nOutputs >= 3) 
    error('Should either have 1 or 2 outputs (just NN or NN and Y).')
end
varargout = cell(1,nOutputs);

if nargin < 4
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

labeler = @(letter,num)sprintf('%s%d_',letter,num);

% Neural network hybrid zonotope
NN = memZono(X,'X');
x0s = arrayfun(@(num){labeler('x_L0_u',num)},1:X.n);
NN.dimKeys = x0s;

prev_xs = x0s;
for i = 1:(length(bs)-1)
    [n2,n1] = size(Ws{i});

    layer = memZono([],[],[],[],[],[]);
    for j = 1:n2
        relu_ij = memZono(relu,sprintf('phi_L%d_u%d_',i,j));
        relu_ij.dimKeys = {sprintf('v_L%d_u%d_',i,j),sprintf('x_L%d_u%d_',i,j)};
        layer = [layer; relu_ij];
    end
    
    vs = arrayfun(@(num){labeler(sprintf('v_L%d_u',i),num)},1:n2);
    xs = arrayfun(@(num){labeler(sprintf('x_L%d_u',i),num)},1:n2);

    % make layer-to-layer connections via constraints
    % v^i = W^{i} x^{i-1} + b*{i}
    V_{i} = NN.affine(bs{i},Ws{i},prev_xs,vs);    
    % NN_lastlayer = Ws{i}*NN.projection(prev_xs);
    % NN_lastlayer.dimKeys = vs; % set the dimKeys so the labeledIntersection finds common dimensions
    % NN_lastlayer = NN_lastlayer + memZono(zono(zeros(length(bs{i}),0),bs{i}),vs);   
    layer = intersect(layer,V_{i},sprintf('intersection_L%i',i));
    % layer = labeledIntersection(layer,NN_lastlayer,vs,sprintf('intersection_L%i',i));
    % %NN = [layer;NN]; % normally this would add many redundant constrains, but memZono identifies exactly repeated constraints
    % NN = layer.cartProd(NN,sprintf('concat_L%i',i));
    % %NN = [NN; labeledIntersection(layer, Ws{i}*NN.projection(prev_xs)+bs{i} ,vs) ];

    NN = NN.intersect(layer,sprintf('concat_L%i',i));

    % Or just:
    % NN = interesect(NN.affine(bs{i},Ws{i},prev_xs,vs), layer);

    

    prev_xs = xs;
end

% make connections from final hidden layer to output
X = NN.projection(x0s);  % this version of X has all the factors from the neural network
Y = Ws{end}*NN.projection(xs);
ys = arrayfun(@(num){labeler('y_',num)},1:length(bs{end}));
Y.dimKeys = ys;
Y = Y + memZono(zono([],bs{end}),ys);

X = X.Z;
Y = Y.Z;

NN = hybZono( ...
                [X.Gc;Y.Gc], ...
                [X.Gb;Y.Gb], ...
                [X.c;Y.c], ...
                X.Ac, ...
                X.Ab, ...
                X.b);

varargout{1} = NN;
if nOutputs == 2
    varargout{2} = Y;
end 

end