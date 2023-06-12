function [x,fVal,exitFlag] = solveMILP(f,A,b,Aeq,beq,lb,ub,vType,opts)

if isempty(opts)
    opts = solverOptions;
end
% Formulated as maximization problem
switch opts.milpSolver
    case 'intlinprog'
%         options = optimoptions(@intlinprog,'Display','off');
%         if opts.nSolutions > 1
%             options.OutputFcn = @savemilpsolutions;
%         end
%         intcon = find(vType == 'B');
%         [x,fVal,exitFlag,optimValues] = intlinprog(-f,intcon,A,b,Aeq,beq,lb,ub,[],options);
    case 'gurobi'
        if ~isempty([A b]) && isempty([Aeq beq])
            model.sense = '<';
            model.A = sparse(A);
            model.rhs = b;
        elseif isempty([A b]) && ~isempty([Aeq beq])
            model.sense = '=';
            model.A = sparse(Aeq);
            model.rhs = beq;
        elseif isempty([A b]) && isempty([Aeq beq])
            model.sense = '=';
            model.A = sparse(zeros(0,length(lb)));
            model.rhs = [];
        else
            error('Can only use < or = constraints (not both) when using gurobi.')
        end
        if ~isempty(f)
            model.obj = f;
        end
        model.lb = lb;
        model.ub = ub;
        model.vtype = vType;
        model.modelsense = 'max';
        params.Threads = 1;
        params.outputflag = 0;
        if opts.nSolutions > 1
            params.PoolSearchMode = 2;
            params.PoolSolutions = opts.nSolutions;
        end
        params.MIPFocus = opts.MIPFocus;
        result = gurobi(model,params);
        if opts.nSolutions == 1
            x = result.x;
            fVal = result.objval;
        else
            if ~isfield(result,'pool')
            	warning('This hybrid zonotope is empty. MILP detected 0 feasible solutions.')
            	x = NaN;
                fVal = NaN;
                exitFlag = -2;
            	return
            end

            x = [result.pool.xn];
            fVal = [result.pool.objval];
        end
        exitFlag = result.status;
    otherwise
        error('Only intlinprog (default) and gurobi are options for solving MILPs.')
end

end