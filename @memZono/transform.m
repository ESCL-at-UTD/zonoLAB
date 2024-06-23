% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       A dimension-aware and memory-preserving transformation with dimensional 
%       specification options.
%   Syntax:
%       out = transform(X,[],M)       <= out = M*X
%       out = transform(X,b,[])       <= out = X + b
%       out = transform(X,b,M)     <= out = M*X + b
%       out = transform(X,Y,[])       <= out = combine(X,Y)
%       out = transform(X,Y,M)     <= out = combine(M*X,Y)
%   Inputs:
%       X       - memZono in R^n
%       b or Y  - vector in R^m or memZono in R^m
%       M       - Transformation matrix R^{m,q}, q<=n
%       inDims  - a cell array of dimensions to be transformed (length q)
%                   if empty, then all dimensions are transformed
%       outDims - a cell array of relabeled dimensions after transformation (length m)
%                   can be empty if matrix is square or a scalar
%   Outputs:
%       Z - memZono in R^m
%           Affine transformation M*X+b is applied to the inDims dimensions
%           If Y memZono is specified, then combine(M*X,Y) is returned
%   Notes:
%       
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function obj = transform(obj1,obj2,M, inDims, outDims)
    arguments
        obj1
        obj2
        M = [];
        inDims = {};
        outDims = {};
    end
    
    %% Input Conditioning

    % if not inDims are specified, then the operation is applied to the entire memZono object.
    if isempty(inDims)
        warning('lack of inDims specification can cause issues with dimension ordering');
        inDims = obj1.dimKeys;
    else
        if ~all(ismember(inDims,obj1.dimKeys)) % checks that if inDims is specified, inDims are all valid labels
            error('inDims specified not all within obj1')
        end
    end

    % outDims must be given unless:
    %   (1) No multiplication --> M is empty
    %   (2) Multiplying, but by a scalar --> M is scalar
    %   (3) Multiplying, but by a square matrix --> M is square
    if isempty(outDims) 
        if ~isempty(M) 
            if isscalar(M) || size(M,1)==size(M,2) % Case 2 or 3
                outDims = inDims;
            else
                error('Output dimension labels (outDims) must be specified for a rectangular matrix multiplication that changes the number of dimensions from input to output.')
            end
        else % Case 1: no multiplication, so dimensions do not change with only the sum
            outDims = inDims;
        end
    % TODO - need to check that if outDims is specified, outDims are all valid labels?
    % (to: @jruths) - this check is dependent on the size of M thus I believe check is done already in other checks (from: @jonaswagner2826)
    end

    % inDims should match the input dimension of the linear map M, if it exists
    if ~isempty(M) && ~isscalar(M) && length(inDims) ~= size(M,2)
        error('The number of inDims labels must match the number of columns of the matrix multiplication.')
    end

    % outDims should match the output dimension of the linear map M, if it exists
    if ~isempty(M) && ~isscalar(M) && length(outDims) ~= size(M,1)
        error('The number of outDims labels must match the number of rows of the matrix multiplication.')
    end

    %% Operation Logic
    if isa(obj2,'memZono') % <-- thus minkowski sum
        obj = combine(affineMap(obj1,[],M,inDims,outDims),obj2); % <-- M Z \oplus Y
    else % <-- thus a vector sum
        obj = affineMap(obj1,obj2,M,inDims,outDims); % <-- M Z + b
    end
end

% Affine Operation - a linear mapping by M and adding a vector b
function out = affineMap(in,b,M,inDims,outDims)
    if isempty(M) && isempty(b)
        out = in;
    elseif isempty(M) % M is empty and the dimensions don't change, so b should match the dimension of the zono already
        out = memZono(in.G, in.c + b, in.A, in.b, in.vset, in.keys);
    else
        if isempty(b)
            b = zeros(length(outDims),1);
        end
        
        % select the dimensions to be mapped by grabbing indices
        % use ismember so that the order of inDims is preserved to match M and b
        [~,M_idx] = ismember(inDims,in.dimKeys);

        % passDims will be the dimensions (keys) that are not going to be mapped
        [~,~,passDims] = memZono.getUniqueKeys(inDims,in.dimKeys);
        % pass_idx are the indices of the passDims keys
        [~,~,~,pass_idx] = memZono.getKeyIndices(inDims,in.dimKeys);

        G_ = [ 
            in.G(pass_idx,:); 
            M*in.G(M_idx,:)
        ];
        c_ = [
            in.c(pass_idx,:);
            M*in.c(M_idx,:) + b
        ];
        % only need to update the dimKeys - reordered and the mapped ones are relabeled
        keys_ = in.keys;
        keys_.dims = [passDims,outDims];
        
        % constraints, constraint keys, factor keys, and vset do not change
        out = memZono(G_, c_, in.A, in.b, in.vset, keys_);
    end
end