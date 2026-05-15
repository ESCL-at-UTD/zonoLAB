function out = reduceInner(obj2,obj1,varargin)
% Method options: 'Hausdorff' or 'Infinity Norm'
disp('Second argument - object that is to be fitted inside the conzono. Can be zono or conzono. Other types PROBABLY not supported')
disp('First argument- Conzono whose inner-approx is desired')


if obj1.n ~= obj2.n
    disp(['Cannot approximate a zonotope in ',num2str(obj2.n),...
        ' dimensions by a zonotope in ',num2str(obj1.n),' dimensions!'])
else
    if isempty(obj1.nC)
        type = 'zono';
    else
        type = 'conZono';
    end
    if isempty(varargin) && strcmp(type,'zono')
        method = 'Hausdorff';
    elseif isempty(varargin) && strcmp(type,'conZono')
        method = 'Infinity Norm';
    else
        method = varargin;
        if  strcmp(method,'Hausdorff') && strcmp(type,'conZono')
            disp('Can only use Hausdorff when inner-approximation is an unconstrained zonotope! Using Infinity Norm instead.')
            method = 'Infinity Norm';
        end
    end
    out = obj1;
    phi_ = sdpvar(obj1.nG,1);
    center_ = sdpvar(obj1.n,1);
    
    obj0 = obj1;
    obj0.c = center_;
    obj0.G = obj1.G*diag(phi_);
    h0 = conZono2AHPoly(obj0);
    h1 = conZono2AHPoly(obj1);
    h2 = conZono2AHPoly(obj2);
    
    cons = [phi_ >= 0];
    opts = sdpsettings('solver','gurobi','verbose',0);
    
    if strcmp(method,'Infinity Norm')
        cons = conContainCheck(h0,h2,cons);
        objs = norm(1e2-phi_,'inf');  %%%%%% Hard coded upper-bound
        
    elseif strcmp(method,'Hausdorff')
        d_ = sdpvar(1,1);
        % B = conZono;
        c = zeros(obj1.n,1);
        G = eye(obj1.n);
        B=zono(G,c);
        B_ah = conZono2AHPoly(B);
        B_d = B_ah*d_;
        
%         h1.c = h0.c;
        h1.c = center_;
        h1.b = diag([phi_;phi_])*h1.b;
        
        cons = [cons, d_ >= 0];
        cons = conContainCheck(h0,h2,cons);
        cons = conContainCheck(h2,h1+B_d,cons);
        objs = d_;  %% Distance alone did not seem to maximize sets.
        objs = objs + norm(1e2-phi_,'inf');  %%%%%% Hard coded upper-bound
    end
    [~] = optimize(cons,objs,opts);
    out.c = value(center_);
    out.G = out.G*diag(value(phi_));
end
end




