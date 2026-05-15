function [R]= bound_r_LP_calculator(G,A,b)
    ng = size(G,2);
    R=zeros(ng,2);    
    for i = 1:ng
        s = zeros(ng,1);
        s(i) = 1;
        lbs=-ones(size(G,2),1);
        lbs(i)=-inf;
        ubs=ones(size(G,2),1);
        ubs(i)=inf;
        lp = Opt('f',-s','Ae',A,'be',b,'lb',lbs,'ub',ubs); % Formulates the linear program
        opt = mpt_solve(lp); % Solves the LP.
        R(i,2) = s'*opt.xopt;
        lp = Opt('f',s','Ae',A,'be',b,'lb',lbs,'ub',ubs); % Formulates the linear program
        opt = mpt_solve(lp); % Solves the LP.
        R(i,1) = s'*opt.xopt;    
    end
end
