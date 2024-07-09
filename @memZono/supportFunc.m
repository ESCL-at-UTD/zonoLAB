% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the value of the support function for a given set and
%       direction
%   Syntax:
%       [s,x] = supportFunc(X,dims,d)
%   Inputs:
%       X - zonotopic set in R^n (memZono object)
%       dims - dimLabels for applying support vector to
%       d - n x 1 support vector direction
%           (n x n_d) support vector array
%   Outputs:
%       s - scalar such that s = max(d'*x), where x \in X
%           (1 x n) array of result value
%       x - n x 1 vector such that x = argmax(d'*x), where x \in X
%           (n x n_d) array of vector results
%   Notes:
%       Defined to find maximum in desired dimension
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [s,x_out] = supportFunc(obj,dims,d_in)
    arguments
        obj memZono
        dims
        d_in = [];
    end

    % Allows selection of keys based on startsWith
    if ~iscell(dims)
        dims = obj.keysStartsWith(dims).dimKeys;
    end

    % Default d_in behavior
    if isempty(d_in)
        d_in = eye(length(dims));
    % (scaler behavior is weird...)
    % elseif isscalar(d_in)
    %     d_in = d_in*eye(length(dims));
    end

    % Construct d
    d = zeros(obj.n,size(d_in,2));
    for i = 1:length(dims)
        idx_dim = strcmp(obj.dimKeys,dims{i});
        for j = 1:size(d_in,2)
            d(idx_dim,j) = d_in(i,j);
        end
    end

    if all([isnumeric(obj),isnumeric(d)],"all")
        % [s,x] = projection(obj,dims).Z.supportFunc(d);
        for j = 1:size(d,2)
            [s(:,j),x(:,j)] = obj.Z.supportFunc(d(:,j));
        end
    else
        if issym(obj)
            error('Sym version not implimented')            
        else %<=== assume optimvar
            % [s,x] = supportFunOptimvar(projection(obj,dims),d);
            % [s,x] = supportFunOptimvar(obj,d);
            for j = 1:size(d,2)
                [s(:,j),x(:,j)] = supportFunOptimvar(obj,d(:,j));
            end
        end
    end

    % Grab only dims of requested result
    for i = 1:length(dims)
        idx_dims(i) = find(strcmp(obj.dimKeys,dims{i}));
    end
    x_out = x(idx_dims,:);


end


function [s,x] = supportFunOptimvar(obj,d)
    switch obj.baseClass
        case 'zono'
            xi = fcn2optimexpr(@(G,d) sign(d'*G)', obj.G, d);
            x = obj.c + obj.G*xi;
            s = d'*x;

        case 'conZono'
            ub = ones(obj.nG,1); lb = -ub;
            xi = fcn2optimexpr(@(G,A,b) ...
                linprog(-(d'*G)',[],[],A,b,lb,ub),...
                obj.G, obj.A, obj.b);
            x = obj.c + obj.G*xi;
            s = d'*x;

        case 'hybZono'
            error('supportFun not currently implimented for optimvar')
        
    end
end