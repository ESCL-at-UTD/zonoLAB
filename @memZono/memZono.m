% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Class:
%       memZono
%       Dimension-aware and memory-encoded Zonotope
%       Z = {c + G \xi | ||\xi||_inf <= 1, A \xi = b, 
%               \xi_\in \{0,1\} \forall_{i \notin vset_}}
%       TODO: finalize definition @jruths - what else do we want for this
%   Syntax:
%       TODO: add
%   Inputs:
%       z - zono object in base zonotope form
%       keys - a struct specifying the underying keys or a string to generate keys from
%       G - n x nG matrix to define zonotope in R^n with nG generators
%       c - n x 1 vector to define center
%       A - nC x nG matrix to define nC equality constraints (A \xi = b)
%       b - nC x 1 vector to define nC equality constraints (A \xi = b)
%       vset_ - nG x 1 logical vector defining if continous or discrete
%   Outputs:
%       Z - memZono object
%   Notes:
%       This class is built upon the functions written for the individual 
%       zonoLAB classes but adds dimensional awareness and memory-encoding/ 
%       preservation within the set operations.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
classdef memZono

    %% Data
    properties (Hidden,Access=private) % Underlying data structure
        G_      % Generator matrix
        c_      % Center
        A_ = [] % Constraint matrix
        b_ = [] % Constraint vector
        vset_    % vset_ defining if generators are continous or discrete
    end

    properties (Dependent,Hidden,Access=private) %(hidden/inaccessble)
        Gc_      % Continuous generator matrix (n x nGc)
        Gb_      % Binary generator matrix (n x nGb)
        Ac_      % Continuous constraint matrix (nC x nGc)
        Ab_      % Binary constraint matrix (nC x nGb)
    end

    properties (Dependent) % These properties get automatically updated when used
        c       % Center (n x 1)
        b       % Constraint vector (nC x 1)
        G       % Generator matrix (n x nG)
        Gc      % Continuous generator matrix (n x nGc)
        Gb      % Binary generator matrix (n x nGb)
        A       % Constraint matrix (nC x nG)
        Ac      % Continuous constraint matrix (nC x nGc)
        Ab      % Binary constraint matrix (nC x nGb)
        vset
    end

    % Dimensions
    properties (Dependent) % These properties get automatically updated when used
        n       % Dimension
        nG      % Number of generators
        nGc     % Number of continuous generators
        nGb     % Number of binary generators
        nC      % Number of constraints
    end
    properties (Dependent,Hidden)
        sz      % [n nG nC]
    end

    % I/O zono
    properties (Dependent,Hidden,Access=private) %<== unable to ensure private... ensure that dimensions are maintained...
        Z_          % Export to a base-zonotope class
    end
    properties (Dependent,Hidden)
        baseClass   % Equivalent base-zonotope class
    end

    % Labeling
    properties (Hidden,Access=protected)
        keys_ = struct( ...
            'factors',[],...
            'dims',[],...
            'cons',[])
    end
    properties (Dependent)
        keys
    end
    properties (Dependent,Hidden)
        factorKeys
        dimKeys
        conKeys
    end

    %% Constructors
    methods
        function obj = memZono(varargin)
            switch nargin
                case 1
                    if isa(varargin{1},'memZono'), obj = varargin{1}; % <--- must have labels
                    else, error('A memZono must be created with labels.')
                    end
                case 2
                    obj.Z_ = varargin{1};
                    obj.keys = varargin{2};
                case 3
                    obj.Z_ = varargin{1};
                    obj.dimKeys = varargin{2};
                    obj.factorKeys = varargin{3};
                    obj.conKeys = varargin{3};
                case 4
                    obj.Z_ = varargin{1};
                    obj.dimKeys = varargin{2};
                    obj.factorKeys = varargin{3};
                    obj.conKeys = varargin{4};
                case 6
                    obj.G_ = varargin{1};
                    obj.c_ = varargin{2};
                    obj.A_ = varargin{3};
                    obj.b_ = varargin{4};
                    obj.vset_ = logical(varargin{5});
                    obj.keys = varargin{6};
                otherwise
                    error('Constructor not specified')
            end
        end
        % Setter function for underlying zonotope data
        function obj = set.Z_(obj,in)
            switch class(in)
                case {'double','sym','optim.problemdef.OptimizationVariable'}
                    obj.c_ = in;
                    obj.vset_ = true(1,0);
                case 'zono'
                    obj.G_ = in.G;
                    obj.c_ = in.c;
                    obj.vset_ = true(1,in.nG);
                case 'conZono'
                    obj.G_ = in.G;
                    obj.c_ = in.c;
                    obj.A_ = in.A;
                    obj.b_ = in.b;
                    obj.vset_ = true(1,in.nG);
                case 'hybZono'
                    obj.G_ = [in.Gc,in.Gb];
                    obj.c_ = in.c;
                    obj.A_ = [in.Ac,in.Ab];
                    obj.b_ = in.b;
                    obj.vset_ = [true(1,in.nGc),false(1,in.nGb)];    
                case 'memZono'
                    obj.G_ = in.G_;
                    obj.c_ = in.c_;
                    obj.A_ = in.A_;
                    obj.b_ = in.b_;
                    obj.vset_ = in.vset_;              
            end
        end
    end
    %% Parameter Set/Read
    methods
        % Matrices
        % Get Matrices  <=== Not For External Use
        function out = get.G_(obj) 
            if isempty(obj.G_); obj.G_ = zeros(obj.n,0); end
            out = obj.G_; 
        end
        function out = get.c_(obj); out = obj.c_; end %<== no checks included
        function out = get.A_(obj)
            if isempty(obj.A_); obj.A_ = zeros(0,obj.nG); end
            out = obj.A_; 
        end
        function out = get.b_(obj); out = obj.b_; end %<== no checks included

        % hybZono Matrices  <=== Not For External Use
        function out = get.Gc_(obj); out = obj.G_(:,obj.vset_); end
        function out = get.Gb_(obj); out = obj.G_(:,~obj.vset_); end
        function out = get.Ac_(obj); out = obj.A_(:,obj.vset_); end
        function out = get.Ab_(obj); out = obj.A_(:,~obj.vset_); end

        % Set hybZono Matrices  <=== Not For External Use
        function obj = set.Gc_(obj,in); obj.G_(:,obj.vset_) = in; end
        function obj = set.Gb_(obj,in); obj.G_(:,~obj.vset_) = in; end
        function obj = set.Ac_(obj,in); obj.A_(:,obj.vset_) = in; end
        function obj = set.Ab_(obj,in); obj.A_(:,~obj.vset_) = in; end

        % Dimensions
        function n = get.n(obj); n = size(obj.c_,1); end
        function nG = get.nG(obj); nG = size(obj.G_,2); end
        function nC = get.nC(obj); nC = size(obj.A_,1); end
        function sz = get.sz(obj); sz = [obj.n,obj.nG,obj.nC]; end
        % hybZono c/b factors
        function nGc = get.nGc(obj); nGc = sum(obj.vset_); end
        function nGb = get.nGb(obj); nGb = sum(~obj.vset_); end
    end
    %% General Methods
    methods
        %% Set Operations (memory-versions)---------------------
        obj = cartProd(obj1,obj2,dims1,dims2,options); % User-facing Cartisian Product
        obj = memoryCartProd(obj1,obj2); % Memory Cartisian Product <=== Not user-facing
        obj = memorySum(obj1,obj2); % Memory Minkowski Sum  (overloaded by plus() )
        obj = memoryIntersection(obj1,obj2,sharedDimLabels); % Memory Intersection Operation (overidden by and())
        % Transformations ---------------------------
        obj = map(obj1,obj2,inDims,outDims); % Maping function
        obj = transform(obj1,obj2,M,inDims,outDims); % Mapping w/ dims (outdated/to be replaced)
        function out = affine(in,M,b,inDims,outDims) %% Affine (technically not needed but useful version instead of map version directly)
            out = in(inDims).transform(b,M,inDims,outDims);
        end
        out = linMap(in,M,inDims,outDims);
        out = funMap(obj1,obj2,inDims,outDims,lbl,funInDims,funOutDims);
        % Plotting -------------------------
        varargout = plot(obj, dims, varargin);
        % Projection & Export ------------------
        out = projection(obj,dims); % projection according to dims
        function out = Z(obj,dims), out = projection(obj,dims).Z_; end
        % Additional methods (implimentations of abstractZono methods in memZono)
        [NN,Y] = reluNN(X,Ws,bs,a);
        function out = boundingBox(in)
            keys_out = in.keys_;
            keys_out.factors = in.dimKeys;
            keys_out.cons = {};
            out = memZono(boundingBox(in.Z_),keys_out);
        end
    end

    %% Indexing  ----------------------------
    methods
        B = subsref(A,S);
        % A = subsasgn(A,S,B); %<---- not overloaded
    end

    %% Overrides ----------------------------
    methods
        % +, plus() - override for memorySum()
        function out = plus(in1,in2), out = in1.memorySum(in2); end
        % *, mtimes() - override for transform (map?) ... mainly just scalers
        function out = mtimes(in1,in2)
            if isa(in2,'memZono'), out = in2.transform([],in1); %<== flip the syntax order
            else, error("mtimes only overloaded for left multiplication");
            end
        end
        % &, and() - overide memoryIntersection w/ checks
        function out = and(obj1,obj2,sharedDimLabels)
            if nargin == 2, error('Must Supply labels for shared dimensions'); end
            out = memoryIntersection(obj1,obj2,sharedDimLabels);
        end
        % Overide for memoryUnion w/ checks ... NOT IMPLIMENTED
        function obj = or(obj1,obj2), error('Union not yet implimented'); end
        % vertcat (Extended cartProd)
        function obj = vertcat(varargin)
            obj = varargin{1};
            for i = 2:nargin, obj = cartProd(obj,varargin{i}); end %<========= Looping approach isn't optimized
        end
        % horzcat not explicity defined... currently returns a cell array of memZono() objects
        function out = horzcat(varargin)
            if nargin == 1, out = varargin{1}; %<= Retun if horzcat not needed
            else, out = varargin(cellfun(@isempty,varargin)); %<== returns as cell array w/ only memZono objects
            end
        end
    end
    methods (Static)
        % sum (extended plus)
        function out = sum(in)
            out = in{1}; for i = 2:numel(in), out = plus(out,in{i}); end %<========= not efficient
        end
        % all (extended and)
        function out = all(in,ptn)
            out = in{1}; if isemtpy(ptn), ptn = 'all'; end
            for i = 2:numel(in), out = and(out,in{i},sprintf('_%s_%d',ptn,i)); end %<========= not efficient
        end
        % any (extended or) ... NOT IMPLIMENTED
        function out = any(in)
            out = in{1}; for i = 2:nargin, out = or(out,in{i}); end
        end
    end
    
    %% Labeling 
    methods
        % Key Getter Functions
        function out = get.keys_(obj); out = obj.keys_; end
        function out = get.keys(obj)
            out = structfun(@(x)(x'), obj.keys_,UniformOutput=false); %<== here for easier visibility in MATLAB editor
        end
        function out = get.dimKeys(obj); out = obj.keys_.dims; end
        function out = get.factorKeys(obj); out = obj.keys_.factors; end
        function out = get.conKeys(obj); out = obj.keys_.cons; end
            
        % Key Setter Functions
        function obj = set.keys(obj,in)
            if isstruct(in)
                obj.keys_ = in;
            else 
                obj.dimKeys = in;
                obj.factorKeys = in;
                obj.conKeys = in;
            end
        end
        function obj = set.dimKeys(obj,in)
            obj.keys_.dims = obj.keysCheck(in,obj.n);
        end
        function obj = set.factorKeys(obj,in)
            obj.keys_.factors = obj.keysCheck(in,obj.nG); 
        end
        function obj = set.conKeys(obj,in)
            obj.keys_.cons = obj.keysCheck(in,obj.nC);
        end
    end
    
    methods
        % Create keys from a patern
        function keys_ = keysStartsWith(obj,pattern)
            for field = string({'dims','factors','cons'})
                keys_.(field) = {};
                for i = 1:length(obj.keys_.(field))
                    if startsWith(obj.keys_.(field){i},pattern)
                        keys_.(field){end+1} = obj.keys_.(field){i};
                    end
                end
            end
        end
        
        % Relabel Dims (takes inDims of obj and returns projected result with outDims)
        function out = relabelDims(obj,inDims,outDims)
            if ~iscell(inDims); inDims = obj.keysStartsWith(inDims).dims; end
            if ~iscell(outDims); outDims = memZono.genKeys(outDims,1:numel(inDims)); end
            out = obj.projection(inDims);
            out.dimKeys = outDims;
        end
    end

    %% Internal Keys Operations 
    methods (Static)
        % Generate keys from prefix  <=== Not For External Use
        function [out] = genKeys(prefix,nums)
            labeler = @(prefix,num)sprintf('%s_%i',prefix,num);
            out = arrayfun(@(num){labeler(prefix,num)},nums);
        end

        % Do Keys Check <=== Not For External Use
        function out = keysCheck(in,n)
            % keysCheck(in,n) - checks to ensure the keys(in) is structured
            % currently for a dimension of n
            if n == 0; out = []; return; end
            if ~iscell(in); in = cellstr(in); end
            if length(in) == n, out = in; 
            elseif isscalar(in)
                out{n} = [];
                for i = 1:n, out{i} = sprintf('%s_%d',in{1},i); end
            elseif length(in)~=n, error('keys not assigned correctly/wrong size');
            elseif numel(unique(out))<numel(out), error('Duplicate keys')
            else, error('keys broken somehow');
            end
        end

        % Get all unique keys  <=== Not For External Use
        function [k1,ks,k2] = getUniqueKeys(in1,in2)
            if isempty(in1) || isempty(in2)
                ks = {}; k1 = in1; k2 = in2;
            else
                ks = intersect(in1,in2);
                k1 = setdiff(in1,ks);
                k2 = setdiff(in2,ks);
            end
        end

        % Finds key indices  <=== Not For External Use
        function [idxk1,idxks1,idxks2,idxk2] = getKeyIndices(in1,in2)
            [k1,ks,k2] = memZono.getUniqueKeys(in1,in2);

            % Find Indices of individual keys
            idxk1 = zeros(1,length(k1));
            idxk2 = zeros(1,length(k2));
            for k = 1:length(k1), idxk1(k) = find(strcmp(in1,k1{k})); end
            for k = 1:length(k2), idxk2(k) = find(strcmp(in2,k2{k})); end

            % Find indices of shared keys
            if isempty(ks)
                idxks1 = []; idxks2 = [];
            else
                idxks1 = zeros(1,length(ks));
                idxks2 = zeros(1,length(ks));
                for k = 1:length(ks)
                    idxks1(k) = find(strcmp(in1,ks{k}));
                    idxks2(k) = find(strcmp(in2,ks{k}));
                end
            end
        end

        % Gets all key labels and indices for dims, factors, and cons  <=== Not For External Use
        function [lbl,idx] = getAllKeyIndices(obj1,obj2)
            %% Keys
            % shared dims
            [lbl.d1,lbl.ds,lbl.d2] = memZono.getUniqueKeys(obj1.dimKeys,obj2.dimKeys);
            [idx.d1,idx.ds1,idx.ds2,idx.d2] = memZono.getKeyIndices(obj1.dimKeys,obj2.dimKeys);
            % shared factors
            [lbl.k1,lbl.ks,lbl.k2] = memZono.getUniqueKeys(obj1.factorKeys,obj2.factorKeys);
            [idx.k1,idx.ks1,idx.ks2,idx.k2] = memZono.getKeyIndices(obj1.factorKeys,obj2.factorKeys);
            % shared cons
            [lbl.c1,lbl.cs,lbl.c2] = memZono.getUniqueKeys(obj1.conKeys,obj2.conKeys);
            [idx.c1,idx.cs1,idx.cs2,idx.c2] = memZono.getKeyIndices(obj1.conKeys,obj2.conKeys);
        end
    end

    %% Input/Output and Display
    methods
        % zono data as tables
        function out = get.G(obj)
            out = array2table(obj.G_, RowNames=obj.dimKeys, VariableNames=obj.factorKeys); 
        end
        function out = get.c(obj)
            out = array2table(obj.c_, RowNames=obj.dimKeys, VariableNames={'c'}); 
        end
        function out = get.A(obj) 
            if obj.nC == 0, out = []; return; end
            out = array2table(obj.A_,RowNames=obj.conKeys,VariableNames=obj.factorKeys); 
        end
        function out = get.b(obj)
            if obj.nC == 0, out = []; return; end
            out = array2table(obj.b_,RowNames=obj.conKeys, VariableNames={'b'}); 
        end
        function out = get.vset(obj)
            out = array2table(reshape(obj.vset_,[],1),RowNames=obj.factorKeys,VariableNames={'vset_'});
        end

        %% Hybzono versions as tables
        function out = get.Gc(obj)
            out = obj.G(:,obj.keys.factors(obj.vset_));
        end
        function out = get.Gb(obj)
            if obj.nGb > 0, out = obj.G(:,obj.keys.factors(~obj.vset_));
            else, out = [];
            end
        end
        function out = get.Ac(obj)
            if obj.nC > 0, out =  obj.A(:,obj.keys.factors(obj.vset_));
            else, out = [];
            end
        end
        function out = get.Ab(obj) 
            if obj.nGb > 0 && obj.nC > 0
                out = obj.A(:,obj.keys.factors(~obj.vset_));
            else, out = [];
            end
        end
        
        % Get base abstractZono class
        function out = get.baseClass(obj)
            if all(obj.vset_)
                if isempty(obj.A_), out = 'zono';
                else, out = 'conZono'; 
                end
            else, out = 'hybZono';
            end
        end

        % Output appropriate base zonotope <=== Not For External Use
        function Z = get.Z_(obj)
            switch obj.baseClass
                case 'zono'
                    Z = zono(obj.G_,obj.c_);
                case 'conZono'
                    Z = conZono(obj.G_,obj.c_,obj.A_,obj.b_);
                case 'hybZono'
                    Z = hybZono(obj.Gc_,obj.Gb_,obj.c_,...
                        obj.Ac_,obj.Ab_,obj.b_);
            end
        end

        % test if special
        function out = isempty(obj)
            out = all([isempty(obj.G_),isempty(obj.c_),...
                isempty(obj.A_),isempty(obj.b_)]);
        end
        function out = issym(obj)
            % tests if any are symbolic
            out = any([isa(obj.G_,'sym'),isa(obj.c_,'sym'),...
                    isa(obj.A_,'sym'),isa(obj.b_,'sym')]);
        end
        function out = isnumeric(obj)
            % tests if all are numeric
            out = all([isnumeric(obj.G_), isnumeric(obj.c_), ...
                    isnumeric(obj.A_), isnumeric(obj.b_)]);
        end
    end
    methods (Static)
        % Returns empty memZono
        function out = empty(); out = memZono([],[],[],[],[],[]); end %<== add dimensions?
    end


    %%%%%%%%%%%%%%%%%%%
    % Questionable Usefulness/Not for Release
    %%%%%%%%%%%%%%%%%%%
    methods (Static)
        % Construction based on name of base objects
        %   (a helper method to combine multiple memZono/abstractZono objects based on name/size)
        function obj = varNameConstructor(varargin)
            for i = 1:nargin
                varName = inputname(i);
                obj_{i} = memZono(varargin{i},varName);
            end
            obj = vertcat(obj_{:});
        end
    end         

    % make these protected/private? --------------------------------------
    methods 
        % Output specific properties
        function varargout = varOut(obj,varargin)
            for i = 1:numel(varargin)
                varargout{i} = obj.(varargin{i});
            end
        end
        % output all data
        function varargout = exportAllData(obj)
            fields = {'G_','c_','A_','b_','vset_','keys_'};
            for i = 1:numel(fields); varargout{i} = obj.(fields{i}); end
        end
    end
end

