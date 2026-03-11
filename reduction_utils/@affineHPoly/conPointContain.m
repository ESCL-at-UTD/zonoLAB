function cons_out = conPointContain(point,obj,cons_in)
% Check if obj1 is contained in obj2
point_n = size(point,1);
if obj.n ~= point_n
    disp(['Cannot check containment of a singleton in ',num2str(point_n),...
        ' dimensions to a AH-polytope in ',num2str(obj.n),' dimensions!'])
else
%     Gamma = sdpvar(obj2.nG,obj1.nG,'full');
%     c_d = 
    beta = sdpvar(obj.nG,1);
    %Lambda =  sdpvar(obj2.nH,obj1.nH,'full');
    
    cons = cons_in;
    %cons = [cons, Lambda >= 0]; % Positive definite
%    cons = [cons, obj1.G == obj2.G*Gamma];
    cons = [cons, obj.c - point == obj.G*beta];
	cons = [cons, abs(beta) <= ones(obj.nG,1)];
%    cons = [cons, Lambda*obj1.A == obj2.A*Gamma];
%    cons = [cons, Lambda*obj1.b <= obj2.b + obj2.A*beta];
    cons_out = cons;
end
end