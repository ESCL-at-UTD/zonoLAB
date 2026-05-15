function [E]= bound_e_LP_calculator(G,A,b)
    ng = size(G,2);
    E=zeros(ng,2);
    
    for i = 1:ng    
        s = zeros(ng,1);
        s(i) = 1;
        lbs=-ones(size(G,2),1);
        ubs=ones(size(G,2),1);
        lp = Opt('f',-s','Ae',A,'be',b,'lb',lbs,'ub',ubs); % Formulates the linear program
        opt = mpt_solve(lp); % Solves the LP.
        E(i,2) = s'*opt.xopt;
        lp = Opt('f',s','Ae',A,'be',b,'lb',lbs,'ub',ubs); % Formulates the linear program
        opt = mpt_solve(lp); % Solves the LP.
        E(i,1) = s'*opt.xopt;    
    end
end
