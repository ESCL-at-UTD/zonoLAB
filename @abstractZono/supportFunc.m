% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the value of the support function for a given set and
%       direction
%   Syntax:
%       [s,x] = supportFunc(X,d)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       d - n x 1 real vector defining direction
%           (N x n) or {N x (n x 1)} can also be passed in to return multiple results at once
%   Outputs:
%       s - scalar such that s = max(d'*x), where x \in X
%       x - n x 1 vector such that x = argmax(d'*x), where x \in X
%   Notes:
%       Defined to find maximum in desired direction
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [s,x] = supportFunc(obj,d)

    obj = obj.simplifySetType;

    if size(d,2) == 1
        x = x_star(obj,d);
        s = d'*x;
    else
        if ~isa(d,'cell'); d_ = num2cell(d',1); else; d_ = d; end
        x_ = cellfun(@(d) x_star(obj,d),d_,"UniformOutput",false);
        s = cellfun(@(d,x) d'*x, d_,x_)';
        if ~isa(d,'cell'); x = horzcat(x_{:})'; else; x = x_; end
        % for i = 1:size(d,1)
        %     x(i,:) = x_star(obj,d(i,:)');
        %     s(i,:) = d(i,:)*x(i,:)';
        % end
    end
end

function x = x_star(obj,d)
    % x_star : computes/solves for the associated x \in X that satisfies max(d'*x)
    
    optSolver = solverOptions; % Using default options

    switch class(obj)
        case 'zono'
            xi = sign(d'*obj.G)';
            x = obj.c + obj.G*xi;
        case 'conZono'
            Aeq = sparse(obj.A);
            beq = [obj.b];
            lb = -ones(obj.nG,1);
            ub =  ones(obj.nG,1);
            [xi,~,~] = solveLP(d'*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
            x = obj.c + obj.G*xi;
        case 'hybZono'
            optSolver.nSolutions = 1;
            optSolver.MIPFocus = 2; % Prove optimality
            % Problem data for mixed-integer linear program (MILP)
            Aeq = [obj.Ac 2*obj.Ab];
            beq = obj.b+obj.Ab*ones(obj.nGb,1);
            lb = -ones(obj.nGc+obj.nGb,1);
            ub =  ones(obj.nGc+obj.nGb,1);
            lb((obj.nGc+1):end) = 0; % Set lower bound on binary factors to zero
            vType(1:obj.nGc) = 'C';
            vType(obj.nGc+1:(obj.nGc+obj.nGb)) = 'B';
            [xi,~,~] = solveMILP(d'*[obj.Gc 2*obj.Gb],[],[],Aeq,beq,lb,ub,vType,optSolver);

            % convert binaries back to {-1 1}
            xi((obj.nGc+1):end,:) = (xi((obj.nGc+1):end,:)-0.5)/0.5;
            x = obj.c + [obj.Gc obj.Gb]*xi;
    end
end