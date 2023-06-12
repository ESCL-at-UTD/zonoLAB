function [v,f] = plotConZono1D(obj,opts)

% Standardized header

Aeq = sparse(obj.A);
beq = [obj.b];
lb = -ones(obj.nG,1);
ub =  ones(obj.nG,1);

dir = 1;
[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,opts);
extreme = [obj.G*x + obj.c]';
v(1,:) = [extreme 0];

dir = -1;
[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,opts);
extreme = [obj.G*x + obj.c]';
v(2,:) = [extreme 0];

f = [1 2];
end