% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example:
%   Reachability of a logical (Boolean) function.
%   Function definition:
%       xi, ui ∈ {0, 1}^20 ,  i ∈ {1,2,3}
%       out1 = u1 ∨ (x2 ⊙ x1)       equiv.  out1 = OR(u1, XNOR(x1, x2))
%       out2 = x2 ⊙ (x1 ∧ u2)       equiv.  out2 = XNOR(x2, AND(x1,u2))
%       out3 = x3,k∼∧ (u2 ⊙ u3)     equiv.  out3 = NAND(x3, XNOR(u2, u3))
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

%% -- Building the hybrid zonotopes for various logical operations
% Firstly, we construct a small little "cube" that we can Minkowski sum to
% the logical sets we want to plot (as those sets originally have only
% discrete points that wouldn't otherwise get plotted properly). 
Add4Plot = hybZono(.1*eye(3),[],[0;0;0],[],[],[]);

% We represent the Boolean logican operations via their graphs, i.e.,
% 3-dimensional sets that represent all pairs of {(input1, input2, output)}

% Construct HZ graph of XNOR
XNOR = hybZono([],[.5;-.5;0],[.5;.5;0],[],[],[]);
XNOR = union(XNOR,hybZono([],[.5;.5;0],[.5;.5;1],[],[],[]));
% plotting version
XNORp = XNOR + Add4Plot;

% Construct HZ graph of AND
AND = hybZono([],[0;-.5;0],[0;.5;0],[],[],[]);
AND = union(AND,hybZono([],[0;.5;.5],[1;.5;.5],[],[],[]));
% plotting version
ANDp = AND + Add4Plot;

% Construct HZ graph of NAND
NAND = diag([1 1 -1])*AND + [0 0 1]';
% plotting version
NANDp = NAND + Add4Plot;

% Construct HZ graph of OR
OR = hybZono([],[0;.5;0],[1;.5;1],[],[],[]);
OR = union(OR,hybZono([],[0;.5;.5],[0;.5;.5],[],[],[]));
% plotting version
ORp = OR + Add4Plot;