function [x,fVal,exitFlag] = solveLP(f,A,b,Aeq,beq,lb,ub,opts)

if isempty(opts)
    opts = solverOptions;
end
% Formulated as maximization problem
switch opts.lpSolver
    case 'linprog'
        options = optimoptions('linprog','Display','off');
        [x,fVal,exitFlag,~] = linprog(-f,A,b,Aeq,beq,lb,ub,options);
    case 'gurobi'
        if ~isempty([A b]) && isempty([Aeq beq])
            model.sense = '<';
            model.A = sparse(A);
            model.rhs = b;
        elseif isempty([A b]) && ~isempty([Aeq beq])
            model.sense = '=';
            model.A = sparse(Aeq);
            model.rhs = beq;
        else
            error('Can only use < or = constraints (not both) when using gurobi.')
        end
        model.obj = f;
        model.lb = lb;
        model.ub = ub;
        model.modelsense = 'max';
        params.Threads = 1;
        params.outputflag = 0;
        result = gurobi(model,params);
        x = result.x;
        fVal = result.objval;
        exitFlag = result.status;
    otherwise
        error('Only linprog (default) and gurobi are options for solving LPs.')
end

end