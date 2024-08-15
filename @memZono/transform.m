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
%       options - Named values to specify additional options
%                 options.retainExtraDims (default = true)  - output includes dims not specified within inDims
%                 options.retainInDims (default = false)  - output includes inDims that are not also specified as outDims
%   Outputs:
%       Z - memZono in R^m
%           Affine transformation M*X+b is applied to the inDims dimensions
%           If Y memZono is specified, then combine(M*X,Y) is returned
%   Notes:
%       
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function obj = transform(obj1,obj2,M, inDims, outDims, options)
    arguments
        obj1
        obj2
        M = [];
        inDims = {};
        outDims = {};
        options.retainExtraDims = true;
        options.retainInDims = false;
    end
    
    %% Input Conditioning

    % if not inDims are specified, then the operation is applied to the entire memZono object.
    if isempty(inDims)
        if ~isempty(M)
            if ~isscalar(M)
                error('lack of inDims specification can cause issues with dimension ordering');
            end
        end
        inDims = obj1.dimKeys;
    else
        if ~all(ismember(inDims,obj1.dimKeys)) % checks that if inDims is specified, inDims are all valid labels
            error('inDims specified not all within obj1')
        end
    end

    % outDims must be given unless:
    %   (1) No multiplication --> M is empty
    %   (2) Multiplying, but by a scalar --> M is scalar
    %   (3) Multiplying, but by a square matrix --> M is square <== we should require this
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

    % if there is no matrix multiplication, then the length of inDims and outDims must be the same
    if isempty(M) && length(inDims) ~= length(outDims)
        error('With no matrix multiplication, the number of inDims and outDims must match.')
    end

    % if a vector sum is done, the length of the vector should match the number of outDims
    if ~isempty(obj2) && ~isa(obj2,'memZono') && length(outDims) ~= size(obj2,1)
        error('The number of outDims must match the length of the vector being added.')
    end

    %% Operation Logic
    if isa(obj2,'memZono') % <-- thus minkowski sum
        obj = combine(affineMap(obj1,[],M,inDims,outDims,options),obj2); % <-- M Z \oplus Y
    else % <-- thus a vector sum
        obj = affineMap(obj1,obj2,M,inDims,outDims,options); % <-- M Z + b
    end
end

% Affine Operation - a linear mapping by M and adding a vector b
function out = affineMap(in,b,M,inDims,outDims,options)
    if isempty(M) % when M is empty, inDims and outDims are equal length
        M = eye(length(inDims));
    end
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
    % if retainExtraDims is false then the pass_idx is eleminated
    if ~options.retainExtraDims; passDims = []; pass_idx = []; end
    % if retainInDims is true then the unspecified inDims
    if options.retainInDims
        retainDims = setdiff(inDims, outDims, "stable");
        [~,~,retain_idx,~]=  memZono.getKeyIndices(retainDims,in.dimKeys);
        passDims = [passDims, retainDims];
        pass_idx = [pass_idx, retain_idx];
    end

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