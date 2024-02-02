function out = affine(obj1,obj2,M, inDims, outDims)

% Syntax:
%   out = Z.affine([],M)       <= out = M Z
%   out = Z.affine(b)       <= out = Z + b
%   out = Z.affine(b,M)     <= out = M Z + b
%   out = Z.affine(Y)       <= out = Z \oplus Y
%   out = Z.affine(Y,M)     <= out = M Z \oplus Y

    arguments
        obj1
        obj2
        M = [];
        inDims = {};
        outDims = {};
    end
    
    %% Input Conditioning
    if isempty(inDims)
        inDims = obj1.dimKeys; 
    end

    %% Operation Logic
    obj1 = obj1.projection(inDims); % select inDims from obj1
    if isempty(obj2) %<== thus linear maping
        out = linMap(obj1,M);
    elseif isa(obj2,'memZono') %<== thus minkowski sum
        if ~isempty(M) %-- transformation matrix
            out = genMinSum(obj1,obj2,M); %<= M Z \oplus Y
        else
            out = minSum(obj1,obj2); %<= Z \oplus Y
        end
    elseif size(obj2,2) == 1 && isempty(M) %<== vector addition
        out = obj1.vecSum(obj2);
    elseif isscalar(M) || size(M,2) == obj1.n  %<== thus Affine Operation
        if isempty(obj2) %<== Linear Map
            out = linMap(obj1,M); %<= M Z
        elseif size(obj2,1) == size(M,1) %<== Affine Operation
            b = obj2; 
            out = affineMap(obj1,M,b); %<= M Z + b
        else
            error("M and b don't line up")
        end
        
    end

    %% Out dimensions
        % @jruths - is this functionality correct? 
        % select specific dims if less then actual?
    if ~isempty(outDims)
        if iscell(outDims) && length(outDims) ~= out.n
            out = out.projection(outDims);
            return; %<== double check correct/desired operation
        end
        out.dimKeys = outDims;
    end




    % if ~isa(obj,'memZono'); error('must be a memZono object'); end % <== not needed?
    % if isempty(M) || isempty(b)
    %     if isempty(M); M = eye(obj.n); end
    %     if isempty(b); b = zeros(size(M,1),1); end
    % elseif isa(b,'memZono')
    
    % Confirm same size
    % TODO - change this check...?  we should just take Minkowski sum of
    % any shared dimensions and the others "pass through"?
    % if obj1.n ~= obj2.n %|| obj1.dimKeys ~= obj2.dimKeys % the second half doesn't work --> ~= is not supported for cell arrays
    %     error('dimensions must be identical')
    % end

    % %% Input Conditioning
    % switch narargin
    %     case 2
    %         if isa(varargin{1},'memZono') %<== minkowski sum
    %             out = minSum(obj,varargin{1});
    %         elseif size(varargin{1},2) == 1 %<== vector sum
    %             out = vecSum(obj,varargin{1});
    %         elseif size(varargin{1},2) == obj.n %<== linear map
    %             out = linMap(obj,varargin{1});
    %         else 
    %             error('wrong inputs/outputs');
    %         end
    %     case 3
    %         if isa(varargin{1},'memZono')
    %             % elseif size(varargin{2},2) == obj.n

    % end





    


    

end

% Linear Map
function out = linMap(in,M)
    out = memZono(M*in.G, M*in.c, in.A, in.b, in.vset, in.keys);
end


% Affine Operation
function out = affineMap(in,M,b)
    out = memZono(M*in.G, M*in.c + b, in.A, in.b, in.vset, in.keys);
end


% Vector Sum
function out = vecSum(in,b)
    out = memZono(in.G, in.c + b, in.A, in.b, in.vset, in.keys);
end

% Generalized Minkowski Sum
function obj = genMinSum(obj1,obj2,M)
    obj = minSum(linMap(obj1,M),obj2);
end

% Minkowski Sum
function obj = minSum(obj1,obj2)
    % get shared keys
    [k1,ks,k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);
    [idxk1,idxks1,idxks2,idxk2] = memZono.getKeyIndices(obj1.factorKeys,obj2.factorKeys);

    % Data selection
    G_ = [obj1.G(:,idxk1),obj1.G(:,idxks1)+obj2.G(:,idxks2), obj2.G(:,idxk2)];
    c_ = obj1.c + obj2.c;
    A_ = [obj1.A(:,idxk1), obj1.A(:,idxks1), zeros(obj1.nC,length(idxk2));
            zeros(obj2.nC,length(idxk1)), obj2.A(:,idxks2), obj2.A(:,idxk2)];
    b_ = [obj1.b; obj2.b];
    
    if obj1.vset(idxks1) ~= obj1.vset(idxks2)
        error('c/d factors not lining up');
        % TODO???: add option for c/d not matching => additonal
        % factors/constraints?
    end
    vset_ = [obj1.vset(idxk1),obj1.vset(idxks1),obj2.vset(idxk2)];

    % Labeling
    keys_.factors = [k1,ks,k2];
    keys_.dims = obj1.dimKeys; %<--- add check for same?
    keys_.cons = [obj1.conKeys; obj2.conKeys];

    obj = memZono(G_,c_,A_,b_,vset_,keys_);
end
