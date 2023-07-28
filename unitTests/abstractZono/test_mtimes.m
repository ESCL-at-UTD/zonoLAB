% Shared variables
% Define set dimensions and matrices
n = 4;
nGc = 15;
nGb = 3;
nC = 8;
Z  = randomSet(1,'zono',n,nGc,nGb,nC);
Zc = randomSet(1,'conZono',n,nGc,nGb,nC);
Zh = randomSet(1,'hybZono',n,nGc,nGb,nC);

rng(1)
M = 1-2*rand(3,4);

% Flag to save new known answers (0 - use old, 1 - save and use new)
saveanswers_mtimes = 0;

if saveanswers_mtimes == 0
    load answers_mtimes.mat
end

%% Test 1: M x zono  
out = M*Z;
if saveanswers_mtimes
    mtimesKnown.out1 = out;
    save answers_mtimes.mat mtimesKnown
end
assert(isequal(out.c,M*Z.c),'Center was not multiplied correctly.')
assert(isequal(out.G,M*Z.G),'Generator matrix was not multiplied correctly.')
assert(isequal(out,mtimesKnown.out1),'Solution does not match known answer.')
%% Test 2: M x conZono  
out = M*Zc;
if saveanswers_mtimes
    mtimesKnown.out2 = out;
    save answers_mtimes.mat mtimesKnown
end
assert(isequal(out.c,M*Zc.c),'Center was not multiplied correctly.')
assert(isequal(out.G,M*Zc.G),'Generator matrix was not multiplied correctly.')
assert(isequal([out.A out.b],[Zc.A Zc.b]),'Constraints were not preserved correctly.')
assert(isequal(out,mtimesKnown.out2),'Solution does not match known answer.')
%% Test 3: M x hybZono  
out = M*Zh;
if saveanswers_mtimes
    mtimesKnown.out3 = out;
    save answers_mtimes.mat mtimesKnown
end
assert(isequal(out.c,M*Zh.c),'Center was not multiplied correctly.')
assert(isequal(out.Gc,M*Zh.Gc),'Continuous generator matrix was not multiplied correctly.')
assert(isequal(out.Gb,M*Zh.Gb),'Binary generator matrix was not multiplied correctly.')
assert(isequal([out.Ac out.Ab out.b],[Zh.Ac Zh.Ab Zh.b]),'Constraints were not preserved correctly.')
assert(isequal(out,mtimesKnown.out3),'Solution does not match known answer.')