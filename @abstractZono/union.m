% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the union of two sets, 
%       Z = X \cup Y = {x \in X || x \in Y}
%   Syntax:
%       Z = union(X,Y)        
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       Y - zonotopic set in R^n (hybZono, conZono, or zono object)
%   Outputs:
%       Z - hybrid zonotope in R^n
%   Notes:
%       Z is always a hybZono.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = union(obj1,obj2)

if obj1.n ~= obj2.n
    error('Inconsistent dimensions.')
end
% Make both input sets hybZono
if ~isa(obj1,'hybZono')
    obj1 = hybZono(obj1);
end
if ~isa(obj2,'hybZono')
    obj2 = hybZono(obj2);
end

% Determine non-zero generators that must be constrained
notZero1c = any(abs(obj1.Gc) > 0);
notZero1b = any(abs(obj1.Gb) > 0);
notZero2c = any(abs(obj2.Gc) > 0);
notZero2b = any(abs(obj2.Gb) > 0);

% Build matrices to enforce or logic within linear equality constraints
% For continuous generators
stair1c = zeros(obj1.nGc,1);	% Initialize at 0
stair1c(notZero1c) = 1;		    % Only constrain non-zeros
stair1c = diag(stair1c);	    % Build matrix
stair1c(~any(stair1c,2),:) = [];% Remove zero rows
stair2c = zeros(obj2.nGc,1);	% Initialize at 0
stair2c(notZero2c) = 1;		    % Only constrain non-zeros
stair2c = diag(stair2c);	    % Build matrix
stair2c(~any(stair2c,2),:) = [];% Remove zero rows

% Top of Ac3
ind1c = sum(notZero1c);
ind2c = sum(notZero2c);
Ac3T = blkdiag([stair1c;-stair1c],[stair2c;-stair2c]);

% For binary generators
stair1b = zeros(obj1.nGb,1);	% Initialize at 0
stair1b(notZero1b) = 1;		    % Only constrain non-zeros
stair1b = diag(stair1b);	    % Build matrix
stair1b( ~any(stair1b,2), : ) = [];	% Remove zero rows
stair2b = zeros(obj2.nGb,1);	% Initialize at 0
stair2b(notZero2b) = 1;		    % Only constrain non-zeros
stair2b = diag(stair2b);	    % Build matrix
stair2b( ~any(stair2b,2), : ) = [];	% Remove zero rows

% Bottom left of Ab3
ind1b = sum(notZero1b);
ind2b = sum(notZero2b);
Ab3B = 0.5*blkdiag([stair1b;-stair1b],[stair2b;-stair2b]);

% Add zeros and 1/2s for binary switch
Ac3 = [ Ac3T ; zeros(2*ind1b+2*ind2b,obj1.nGc+obj2.nGc) ];
Ab3 = [ zeros(2*ind1c+2*ind2c,obj1.nGb+obj2.nGb) ; Ab3B ];
Ab3 = [ Ab3 , [ 0.5*ones(2*ind1c,1) ; -0.5*ones(2*ind2c,1) ; 0.5*ones(2*ind1b,1) ; -0.5*ones(2*ind2b,1) ] ];

% The b vector
b3 = [ 0.5*ones(2*ind1c,1) ; 0.5*ones(2*ind2c,1) ; zeros(ind1b,1) ; ones(ind1b,1) ; zeros(ind2b,1) ; ones(ind2b,1) ];

% Find hat variables for substitutions 
if obj2.nGb == 0
	c1 = obj1.c;
else
	c1 = obj2.Gb*ones(obj2.nGb,1)+obj1.c;
end
if obj1.nGb == 0
	c2 = obj2.c;
else
	c2 = obj1.Gb*ones(obj1.nGb,1)+obj2.c;
end
G_hat = [ 0.5*eye(obj1.n) -0.5*eye(obj1.n) ; 0.5*eye(obj1.n) 0.5*eye(obj1.n) ]*[c1;c2];
Gb_hat = G_hat(1:obj1.n,:);
c_hat = G_hat(obj1.n+1:end,:);

A1_hat = [ -0.5*eye(obj1.nC) 0.5*eye(obj1.nC) ; 0.5*eye(obj1.nC) 0.5*eye(obj1.nC) ]*[obj1.b;-1*obj1.Ab*ones(obj1.nGb,1)];
A1b_hat = A1_hat(1:obj1.nC,:);
b1_hat = A1_hat(obj1.nC+1:end,:);

A2_hat = [ -0.5*eye(obj2.nC) 0.5*eye(obj2.nC) ; 0.5*eye(obj2.nC) 0.5*eye(obj2.nC) ]*[-1*obj2.Ab*ones(obj2.nGb,1);obj2.b];
A2b_hat = A2_hat(1:obj2.nC,:);
b2_hat = A2_hat(obj2.nC+1:end,:);

% Put it all together
Gc = [ obj1.Gc obj2.Gc zeros(obj1.n,2*ind1c+2*ind1b+2*ind2c+2*ind2b) ];
Gb = [ obj1.Gb obj2.Gb Gb_hat ];
c = c_hat;
Ac = blkdiag(obj1.Ac,obj2.Ac);
Ac = [ Ac zeros(obj1.nC+obj2.nC,(2*ind1c+2*ind1b+2*ind2c+2*ind2b)); Ac3 eye(2*ind1c+2*ind1b+2*ind2c+2*ind2b) ];
Ab = [ obj1.Ab zeros(obj1.nC,obj2.nGb) A1b_hat ; zeros(obj2.nC,obj1.nGb) obj2.Ab A2b_hat ];
Ab = [ Ab ; Ab3 ];
b = [ b1_hat ; b2_hat ; b3 ];
out = hybZono(Gc,Gb,c,Ac,Ab,b);

end