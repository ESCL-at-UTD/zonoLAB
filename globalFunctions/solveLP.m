% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Global function:
%       Returns the maximum value x \in R^n of a linear program (LP)
%   Syntax:
%       [x,fVal,exitFlag] = solveLP(f,A,b,Aeq,beq,lb,ub,opts)
%   Inputs:
%       f - n x 1 vector defining linear objective function to be maximized
%       A - nC x n matrix defining equality constraints (A x = b)
%       b - nC x 1 vector defining eqaulity constraints (A x = b)
%       Aeq - nCeq x n matrix defining inequality constraints (Aeq x <= beq)
%       beq - nCeq x 1 vector defining ineqaulity constraints (Aeq x <= beq)
%       lb - n x 1 vector defining lower bounds (lb <= x)
%       ub - n x 1 vector defining upper bounds (x <= ub)
%       optSolver - solver options needed for linear propgram
%   Outputs:
%       x - n x 1 vector maximizing objective function subject to constraints
%       fVal - 1 x 1 scalar maximum value of objective function
%       exitFlag - 1 x 1 scalar exit condition (see solver-specific exit condition codes)
%   Notes:
%       Currently, only linprog and Gurobi are supported.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [x,fVal,exitFlag] = solveLP(f,A,b,Aeq,beq,lb,ub,optSolver)

if isempty(optSolver)
    optSolver = solverOptions;
end
% Formulated as maximization problem
switch optSolver.lpSolver
    case 'linprog'
        options = optimoptions('linprog','Display','off');
        [x,fVal,exitFlag,~] = linprog(-f,A,b,Aeq,beq,lb,ub,options);
    case 'gurobi'
        if ~isempty([A b]) && isempty([Aeq beq])
            model.sense = '<';
            model.A = sparse(A);
            model.rhs = double(full(b));
        elseif isempty([A b]) && ~isempty([Aeq beq])
            model.sense = '=';
            model.A = sparse(Aeq);
            model.rhs = double(full(beq));
        elseif ~isempty([A b]) && ~isempty([Aeq beq]) % Need to fix case where both are empty
            error('Can only use < or = constraints (not both) when using gurobi.')
        end
        model.obj = f;
        model.lb = lb;
        model.ub = ub;
        model.modelsense = 'max';
        params.Threads = 1;
        params.outputflag = 0;
        result = gurobi(model,params);
        exitFlag = result.status;
        x = NaN;
        fVal = NaN;
        if strcmp(exitFlag,'OPTIMAL')
            x = result.x;
            fVal = result.objval;
        else
            warning('Linear program did not find an optimal solution.');
        end
    otherwise
        error('Only linprog (default) and gurobi are options for solving LPs.')
end

end