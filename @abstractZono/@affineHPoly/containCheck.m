function out = containCheck(obj1,obj2)
% Check if obj1 is contained in obj2
if obj1.n ~= obj2.n
    disp(['Cannot check containment of a AH-polytope in ',num2str(obj1.n),...
        ' dimensions to a AH-polytope in ',num2str(obj2.n),' dimensions!'])
else
    Gamma = sdpvar(obj2.nG,obj1.nG,'full');
    beta = sdpvar(obj2.nG,1);
    Lambda =  sdpvar(obj2.nH,obj1.nH,'full');
    
    cons = [];
    cons = [cons, Lambda >= 0]; % Positive definite
    cons = [cons, obj1.G == obj2.G*Gamma];
    cons = [cons, obj2.c - obj1.c == obj2.G*beta];
    cons = [cons, Lambda*obj1.A == obj2.A*Gamma];
    cons = [cons, Lambda*obj1.b <= obj2.b + obj2.A*beta];
    
    objs = 0; % Feasibility problem
    opts = sdpsettings('solver','linprog','verbose',0); % Solver settings
    sol = optimize(cons,objs,opts); % Solves the LP.
    if sol.problem == 0
        out = 1;
    else
        out = 0;
    end
end
end