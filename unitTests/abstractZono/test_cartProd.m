% Shared variables
% Define set 1 dimensions and matrices
n1 = 2;
nGc1 = 15;
nGb1 = 3;
nC1 = 8;
Z1  = randomSet(1,'zono',n1,nGc1,nGb1,nC1);
Zc1 = randomSet(1,'conZono',n1,nGc1,nGb1,nC1);
Zh1 = randomSet(1,'hybZono',n1,nGc1,nGb1,nC1);
% Define set 2 dimensions and matrices
n2 = 3;
nGc2 = 10;
nGb2 = 4;
nC2 = 5;
Z2  = randomSet(2,'zono',n2,nGc2,nGb2,nC2);
Zc2 = randomSet(2,'conZono',n2,nGc2,nGb2,nC2);
Zh2 = randomSet(2,'hybZono',n2,nGc2,nGb2,nC2);

% Flag to save new known answers (0 - use old, 1 - save and use new)
saveanswers_cartProd = 0;

if saveanswers_cartProd == 0
    load answers_cartProd.mat
end

%% Test 1: zono x zono  
out = cartProd(Z1,Z2);
if saveanswers_cartProd
    cartProdKnown.out1 = out;
    save answers_cartProd.mat cartProdKnown
end
assert(isequal(out,cartProdKnown.out1),'Solution does not match known answer.')

%% Test 2: conZono x conZono  
out = cartProd(Zc1,Zc2);
if saveanswers_cartProd
    cartProdKnown.out2 = out;
%     save answers_cartProd.mat cartProdKnown
end
assert(isequal(out,cartProdKnown.out2),'Solution does not match known answer.')
%% Test 3: hybZono x hybZono  
out = cartProd(Zh1,Zh2);
if saveanswers_cartProd
    cartProdKnown.out3 = out;
    save answers_cartProd.mat cartProdKnown
end
assert(isequal(out,cartProdKnown.out3),'Solution does not match known answer.')
%% Test 4: zono x conZono  
out = cartProd(Z1,Zc2);
if saveanswers_cartProd
    cartProdKnown.out4 = out;
    save answers_cartProd.mat cartProdKnown
end
assert(isequal(out,cartProdKnown.out4),'Solution does not match known answer.')
%% Test 5: conZono x zono  
out = cartProd(Zc1,Z2);
if saveanswers_cartProd
    cartProdKnown.out5 = out;
    save answers_cartProd.mat cartProdKnown
end
assert(isequal(out,cartProdKnown.out5),'Solution does not match known answer.')
%% Test 6: zono x hybZono  
out = cartProd(Z1,Zh2);
if saveanswers_cartProd
    cartProdKnown.out6 = out;
    save answers_cartProd.mat cartProdKnown
end
assert(isequal(out,cartProdKnown.out6),'Solution does not match known answer.')
%% Test 7: hybZono x zono  
out = cartProd(Zh1,Z2);
if saveanswers_cartProd
    cartProdKnown.out7 = out;
    save answers_cartProd.mat cartProdKnown
end
assert(isequal(out,cartProdKnown.out7),'Solution does not match known answer.')
%% Test 8: conZono x hybZono  
out = cartProd(Zc1,Zh2);
if saveanswers_cartProd
    cartProdKnown.out8 = out;
    save answers_cartProd.mat cartProdKnown
end
assert(isequal(out,cartProdKnown.out8),'Solution does not match known answer.')
%% Test 9: hybZono x conZono  
out = cartProd(Zh1,Zc2);
if saveanswers_cartProd
    cartProdKnown.out9 = out;
    save answers_cartProd.mat cartProdKnown
end
assert(isequal(out,cartProdKnown.out9),'Solution does not match known answer.')