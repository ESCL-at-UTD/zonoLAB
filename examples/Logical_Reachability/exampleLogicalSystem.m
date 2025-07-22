% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
% Example:
%   Reachability of a logical (Boolean) function.
%   Function definition:
%       xi, ui ∈ {0, 1}^20 ,  i ∈ {1,2,3}
%       out1 = u1 ∨ (x2 ⊙ x1)       equiv.  out1 = OR(u1, XNOR(x1, x2))
%       out2 = x2 ⊙ (x1 ∧ u2)       equiv.  out2 = XNOR(x2, AND(x1,u2))
%       out3 = x3 ∼∧ (u2 ⊙ u3)     equiv.  out3 = NAND(x3, XNOR(u2, u3))
%
% For more details, see the following papers:
% 
% A. Alanwar et al., "Polynomial Logical zonotope: A Set Representation 
% for Reachability Analysis of Logical Systems," in Automatica, 171, 2025.
%
% J. A. Siefert et al., "Reachability Analysis Using Hybrid Zonotopes 
% and Functional Decomposition," in IEEE Transactions on Automatic Control
% (appearing soon).
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

%% -- Constructing the state-update set for the function of interest
% This is the function from "Polynomial Logical Zonotopes: A Set
% Representation for Reachability Analysis of Logical Systems" by Alanwar
% et al., and was also examined in the recent TAC paper by Siefert et al.

% We construct a HZ representing the graph of the function by enforcing the
% sets constructed earlier along the corresponding dimensions of the basis set.

% Create basis
Phi = hybZono(zeros(12,1),.5*eye(12),.5*ones(12,1),[],[],[]);
% Couple the intermediate variables
% [b1, b2, b3, u1, u2, u3, tmp1=xnor(b2,b1), tmp2=and(b1,u2), tmp3=xnor(u2,u3), b1+, b2+, b3+]^T
Phi = and(Phi,XNOR,M_i([2 1 7],12));
Phi = and(Phi,AND ,M_i([1 5 8],12));
Phi = and(Phi,XNOR,M_i([5 6 9],12));
% Couple k+1 variables
% [b1, b2, b3, u1, u2, u3, tmp1, tmp2, tmp3, b1+=or(u1,tmp1), b2+=xnor(b2,tmp2), b3+=nand(b3,tmp3)]^T
Phi = and(Phi,  OR,M_i([4 7 10],12));
Phi = and(Phi,XNOR,M_i([2 8 11],12));
Phi = and(Phi,NAND,M_i([3 9 12],12));
% Tranform to get rid of intermediate states
Phi = M_i([1:6, 10:12],12)*Phi;

% complexity reduction
Phi = removeRedundancy(Phi);

%% -- Testing the state-update set
% We test the state-update set by first randomly generating points in the 
% input space and representing each as a HZ. We intersect each point with
% the HZ state-update set to obtain a single-point output set. To ensure we
% constructed Phi correctly, we compare that point with the expected output
% of the function, and find that it works correctly in all tested cases. 
errors_found = 0;
num_test_cases = 100;
clear Phi_test_data
faulty_cases = [];
for i = 1:num_test_cases
    % generate random "logical" points and represent as HZ
    R3{i} = round(rand(3,1));
    R3set = hybZono([],[],R3{i},[],[],[]);
    U{i} = round(rand(3,1));
    Uset = hybZono([],[],U{i},[],[],[]);
    
    % Successor Set Identity
    R3plusSet = [zeros(3,6) eye(3)]*and(Phi,cartProd(R3set,Uset),[eye(6) zeros(6,3)]);
    
    % Note Gc=0, and # leaves = 1
    R3plus = R3plusSet.c+R3plusSet.Gb*R3plusSet.getLeaves;
    % convert to MATLAB logical variables
    R3{i} = logical(R3{i});
    U{i} = logical(U{i});
    
    % Should anything go wrong, you can uncover the faulty cases here
    Phi_test_data{i} = [R3{i}, U{i}, R3plus, BooleanFunction(R3{i},U{i})];
    
    temp = max(abs(R3plus-BooleanFunction(R3{i},U{i})));
    if temp >= 1e-10 % rounding errors happen
        faulty_cases = [faulty_cases i];
        errors_found = errors_found + max(abs(R3plus-BooleanFunction(R3{i},U{i})));
    end
end
errors_found

%% -- Stack the system to make it high-dimensional
% In the TAC paper, we stack this system upon itself 20 times to test our
% methods in high-dimensions. 

Phi20 = Phi;
for i = 2:20
   Phi20 = cartProd(Phi20,Phi); 
end

% Construct linear tranform to untangle states, inputs, and updated states
order = [];
for i = 1:9
    order = [order i:9:180];
end
T = M_i(order,180);

% Untangle states, inputs, and updated states
Phi20 = T*Phi20;

%% -- Performing reachability analysis using the state-update set
% Parameters for the reachability analysis:
num_steps = 5;
numOfVar = 3;
numOfVarU = 3;
numOfPoints = 8;

% Generate random points in the X input space
Apt = randi([0 1], numOfVar*20,  numOfPoints);
% Initial set represented as a hybrid zonotope
R{1} = hybZono({Apt,eye(numOfPoints)});

% Generate random points in the U space
Upt = randi([0 1], numOfVarU*20, numOfPoints);
U = hybZono({Upt, eye(numOfPoints)});
% All possible points in the U input space
% U = hybZono(zeros(numOfVarU*20,1),diag(.5*ones(60,1)),.5*ones(60,1),[],[],[]);

ZeroEye = [zeros(60,120)    eye(60)];
EyeZero = [eye(120)         zeros(120,60)];
clear N_StepTime_TotalTime
tic
for i = 1:num_steps
   startIt = toc;
   % Successor Set Identity
   R{i+1} = ZeroEye*and(Phi20,cartProd(R{i},U),EyeZero);
   % Keeping track of computation times
   N_StepTime_TotalTime(i,:) = [i toc-startIt, toc];
end

% Display computation times
N_StepTime_TotalTime


%% -- Comparison to Logical Zonotopes
% Since the writing of the TAC paper, some improvements have been made to
% the toolboxes for both hybrid zonotopes and logical zonotopes (which can
% be found here: https://github.com/aalanwar/Logical-Zonotope.git). So, now
% we code this up to utilize the most recent capabilities of each toolbox. 

% To run this portion of the code, please ensure the folder
% 'Logical-Zonotope' is added to the path. 

% The following portion of code is copied from the file
% 'Logical-Zonotope/BoolFunctionExample.m'.

clear

steps_LZ = 5; % number of steps in the prediction horizon
dim_stack = 20; % number of times we "stacked" the system on itself
numOfVar = 3;
numOfVarU = 3;

numOfGens = 4;
numOfGensU =4;

for i=1:numOfVar 
    c= logical(randi([0 1],dim_stack,1));
    g ={};
    for j=1:numOfGens
        g{j} = logical(randi([0 1],dim_stack,1));
    end
    Az{i} = logicalZonotope(c,g);
    Apz{i} = logicalPolyZonotope(c,g,eye(length(g)));
    Apt{i} = evaluate(Az{i});
end

for i=1:numOfVarU 
    c= logical(randi([0 1],dim_stack,1) );
    g ={};
    for j=1:numOfGensU
        g{j} = logical(randi([0 1],dim_stack,1));
    end
    Uz{i} = logicalZonotope(c,g);
    Upz{i} = logicalPolyZonotope(c,g,eye(length(g)));
    pointsU{i} = evaluate(Uz{i});
end

tic
for i =1:steps_LZ
    [Az] = customBoolFun(Az,Uz);
end
execlogicZono = toc

tic
for i =1:steps_LZ
    [Apz] = customBoolFun(Apz,Upz);
end
execplogicPolyZono = toc

%% ------------------SUBFUNCTIONS-------------------------

function [e] = e_i(i,n)
    e = zeros(1,n);
    e(i) = 1;
end

function [M] = M_i(inds,n)
    M = [];
    for i = 1:length(inds)
        e = e_i(inds(i),n);
        M = [M;e];
    end
end

function out = BooleanFunction(in,u)
    out = [ or(u(1),~xor(in(2),in(1)));
            ~xor(in(2),and(in(1),u(2)));
            ~and(in(3),~xor(u(2),u(3)))];
end

function A = customBoolFun(A,U)
temp{1} = or(U{1},xnor(A{2},A{1}));
temp{2} = xnor(A{2},and(A{1},U{2}));
xn = xnor(U{2},U{3});
temp{3} = nand(A{3},xn);
A = temp;
end