% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return vertices and faces for a constrained zonotope in 1D
%   Syntax:
%       [v,f] = plotConZono1D(Z,optSolver)
%   Inputs:
%       Z - 1D constrained zonotope in CG-Rep (conZono object)
%       optSolver - solver options needed for linear propgram
%   Outputs:
%       v - 2 x 2 matrix, first row indicating lower and upper bounds (x-direction),
%                         second row indicting zeros for plotting purposes (y-direction)
%       f - 1 x 2 vector, indicating a single face containing both vertices 
%   Notes:
%       Not intended to be called directly by user.
%       Use [v,f] = plot(obj,varargin) instead (method of abstractZono)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [v,f] = plotConZono1D(obj,optSolver)

% Problem data for linear program (LP)
Aeq = sparse(obj.A);
beq = [obj.b];
lb = -ones(obj.nG,1);
ub =  ones(obj.nG,1);

try
    dir = 1; % Find upper bound
    %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
    x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
    extreme = [obj.G*x + obj.c]';
    v(1,:) = [extreme 0];

    dir = -1; % Find lower bound
    %[x,~,~] = solveLP(dir*obj.G,[],[],Aeq,beq,lb,ub,optSolver);
    x = findVertex(dir,obj.G,Aeq,beq,lb,ub,optSolver);
    extreme = [obj.G*x + obj.c]';
    v(2,:) = [extreme 0];

    f = [1 2];
catch E
    if strcmp(E.identifier, 'PlotError:VertexNotFound')
        if checkEmpty(conZono(obj.G,obj.c,Aeq,beq))
            warning('zonoLAB:EmptyZonotope','Constrained zonotope is empty and cannot be plotted.')
            v = []; f = [];
            return
        end
        rethrow(E);
    else
        rethrow(E);
    end
end

end

% Local functions
function [x] = findVertex(dir,G,Aeq,beq,lb,ub,optSolver)
    [x,~,~] = solveLP(dir*G,[],[],Aeq,beq,lb,ub,optSolver);
    if isnan(x)
        error('PlotError:VertexNotFound','Could not find a solution for a vertex while plotting')
    end
end