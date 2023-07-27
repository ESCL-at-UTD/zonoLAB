% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Global function:
%       Returns the maximum value x \in R^n of a mixed-integer linear program (MILP)
%   Syntax:
%       [x,fVal,exitFlag] = solveMILP(f,A,b,Aeq,beq,lb,ub,vType,opts)
%   Inputs:
%       f - n x 1 vector defining linear objective function to be maximized
%       A - nC x n matrix defining equality constraints (A x = b)
%       b - nC x 1 vector defining eqaulity constraints (A x = b)
%       Aeq - nCeq x n matrix defining inequality constraints (Aeq x <= beq)
%       beq - nCeq x 1 vector defining ineqaulity constraints (Aeq x <= beq)
%       lb - n x 1 vector defining lower bounds (lb <= x)
%       ub - n x 1 vector defining upper bounds (x <= ub)
%       vType - n x 1 character vector denoting if each element of x is a
%               continuous ('C') or binary ('B') variable
%       optSolver - solver options needed for mixed-integer linear propgram
%   Outputs:
%       x - n x 1 vector maximizing objective function subject to constraints
%       fVal - 1 x 1 scalar maximum value of objective function
%       exitFlag - 1 x 1 scalar exit condition (see solver-specific exit condition codes)
%   Notes:
%       Currently, only Gurobi is supported.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [x,fVal,exitFlag] = solveMILP(f,A,b,Aeq,beq,lb,ub,vType,optSolver)

if isempty(optSolver)
    optSolver = solverOptions;
end
% Formulated as maximization problem
switch optSolver.milpSolver
    case 'intlinprog'
        error('Only gurobi is supported at this time.')
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
        if optSolver.nSolutions > 1
            params.PoolSearchMode = 2;
            params.PoolSolutions = min(optSolver.nSolutions,2e9);
        end
        params.MIPFocus = optSolver.MIPFocus;
        result = gurobi(model,params);
        if optSolver.nSolutions == 1
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