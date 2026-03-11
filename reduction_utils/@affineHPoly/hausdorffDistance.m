function out = hausdorffDistance(obj1,obj2)
if obj1.n ~= obj2.n
    disp(['Cannot compute the Hausdorff distance of a AH-polytope in ',num2str(obj1.n),...
        ' dimensions to a AH-polytope in ',num2str(obj2.n),' dimensions!'])
else
    d_ = sdpvar(1,1);
    B = conZono;
    B.c = zeros(obj1.n,1);
    B.G = eye(obj1.n);
    B_ah = conZono2AHPoly(B);
    B_d = B_ah*d_;
    
    cons = [d_ >= 0];
    cons = conContainCheck(obj1,obj2+B_d,cons);
    cons = conContainCheck(obj2,obj1+B_d,cons);
    
    objs = d_; % minimize Hausdorff distance
    opts = sdpsettings('solver','linprog','verbose',0); % Solver settings
    sol = optimize(cons,objs,opts); %#ok<NASGU> % Solves the LP.
    out = value(d_);
end
end