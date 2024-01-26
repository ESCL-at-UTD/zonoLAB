% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Class:
%       memZono
%       Memory-based Zonotope
%   Syntax:
%       TODO: add
%   Inputs:
%       TODO: add
%   Outputs:
%       TODO: add
%   Notes:
%       This class is built upon the functions written for the individual 
%       zonoLAB classes but adds memory functionality and associated set 
%       operations.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
classdef memZono

    %% Data
    properties (Hidden)
        G_
        c_
        A_ = []
        b_ = []
        vset
        % tags struct = struct( ...
        %     'factors',struct([]),...
        %     'dims',struct([]),...
        %     'cons',struct([]))
    end

    properties (Dependent,Hidden)
        G
        A
    end

    properties (Dependent)
        c
        b
        % G
        Gc
        Gb
        % A
        Ac
        Ab
    end

    % I/O zono
    properties (Dependent, Hidden)
        Z
        baseClass
    end

    % Dimensions
    properties (Dependent)
        n
        nG
        nGc
        nGb
        nC
    end

    % Labeling
    properties (Hidden)
        keys = struct( ...
            'factors',[],...
            'dims',[],...
            'cons',[])
    end

    properties (Dependent) % currently not hidden and simplified to just keys for now
        factorKeys
        dimKeys
        conKeys
    end

    %% Constructor
    methods
        function obj = memZono(varargin)
            if nargin == 2
                obj.Z = varargin{1};
                obj.keys = varargin{2};

                % if isstruct(varargin{2})
                %     obj.tags = varargin{2};
                % elseif isalpha_num(varargin{2})
                %     obj.factorTags = varargin{2};
                % end

                % obj.tags = varargin{2};
                % switch class(varargin{1})
                %     case 'hybZono'
                % end
                
                % obj.factorKeys = varargin{2};
                % obj.dimKeys = varargin{2};
                % obj.dimKeys = varargin{2};
            elseif nargin == 6
                obj.Z = conZono(varargin{1:4});
                % obj.G_ = varargin{1};
                % obj.c_ = varargin{2};
                % obj.A_ = varargin{3};
                % obj.b_ = varargin{4};
                %  TODO: add checks?
                obj.vset = varargin{5};
                % obj.factorKeys = varargin{end};
                obj.keys = varargin{6};
            else
                error('non-simple constructor not specified')
            end
        end
    end

    %% Get/Set Functions
    methods
        %% System Definitions
        
        % Standard Matrices
        function G = get.G(obj); G = obj.G_; end
        function c = get.c(obj); c = obj.c_; end
        function A = get.A(obj)
            if isempty(obj.A_); obj.A_ = zeros(0,obj.nG); end
            A = obj.A_; 
        end
        function b = get.b(obj); b = obj.b_; end

        function obj = set.G(obj,G); obj.G_ = G; end
        function obj = set.c(obj,c); obj.c_ = c; end
        function obj = set.A(obj,A); obj.A_ = A; end
        function obj = set.b(obj,b); obj.b_ = b; end

        % HybZono Matrices
        function Gc = get.Gc(obj); Gc = obj.G(:,obj.vset); end
        function Gb = get.Gb(obj); Gb = obj.G(:,~obj.vset); end
        function Ac = get.Ac(obj); Ac = obj.A(:,obj.vset); end
        function Ab = get.Ab(obj); Ab = obj.A(:,~obj.vset); end

        % TODO: set for hyb stuff...

        % Dimensions
        function n = get.n(obj); n = size(obj.G,1); end
        function nG = get.nG(obj); nG = size(obj.G,2); end
        function nGc = get.nGc(obj); nGc = size(obj.Gc,2); end
        function nGb = get.nGb(obj); nGb = size(obj.Gb,2); end
        function nC = get.nC(obj); nC = size(obj.A,1); end
        
        % In/Out Zono
        function out = get.baseClass(obj)
            if all(obj.vset)
                if isempty(obj.A_)
                    out = 'zono';
                else
                    out = 'conZono';
                end
            else
                out = 'hybZono';
            end
        end

        function Z = get.Z(obj)
            switch obj.baseClass
                case 'zono'
                    Z = zono(obj.G,obj.c);
                case 'conZono'
                    Z = conZono(obj.G,obj.c,obj.A,obj.b);
                case 'hybZono'
                    Z = hybZono(obj.Gc,obj.Gb,obj.c,...
                        obj.Ac,obj.Ab,obj.b);
            end
        end

        function obj = set.Z(obj,Z)
            switch class(Z)
                case 'zono'
                    obj.G_ = Z.G;
                    obj.c_ = Z.c;
                    obj.vset = true(1,Z.nG);
                case 'conZono'
                    obj.G_ = Z.G;
                    obj.c_ = Z.c;
                    obj.A_ = Z.A;
                    obj.b_ = Z.b;
                    obj.vset = true(1,Z.nG);
                case 'hybZono'
                    obj.G_ = [Z.Gc,Z.Gb];
                    obj.c_ = Z.c;
                    obj.A_ = [Z.Ac,Z.Ab];
                    obj.b_ = Z.b;
                    obj.vset = [true(1,Z.nGc),false(1,Z.nGb)];
            end
        end

    end
    
    %% Labeling
    methods
        function out = get.keys(obj); out = obj.keys; end
        function out = get.factorKeys(obj); out = obj.keys.factors; end
        function out = get.dimKeys(obj); out = obj.keys.dims; end
        function out = get.conKeys(obj); out = obj.keys.cons; end

        function obj = set.keys(obj,in)
            if isstruct(in) %<-- add better checks?
                obj.keys = in;
            else
                obj.factorKeys = in;
                obj.dimKeys = in;
                obj.conKeys = in;
            end
        end
        function obj = set.factorKeys(obj,in)
            try obj.keys.factors = obj.keysCheck(in,obj.nG); 
            catch; warning('factor key set issue');
            end
        end
        function obj = set.dimKeys(obj,in)
            try obj.keys.dims = obj.keysCheck(in,obj.n);  
            catch; warning('dim key set issue'); 
            end
        end
        function obj = set.conKeys(obj,in)
            try obj.keys.cons = obj.keysCheck(in,obj.nC); 
            catch; warning('con key set issue');
            end
        end
           
    end

    methods (Static)
        function out = keysCheck(in,n)
            if n == 0; out = []; return; end
            if ~iscell(in); in = cellstr(in); end
            if length(in) == n
                out = in; 
            elseif length(in) == 1
                out{n} = [];
                for i = 1:n
                    out{i} = sprintf('%s_%d',in{1},i);
                end
            elseif length(in) ~= n
                warning('keys not assigned correctly/wrong size')
            else
                error('keys broken');
            end
        end

        function [k1,ks,k2] = getUniqueKeys(in1,in2)
            if isempty(in1) || isempty(in2)
                ks = {};
                k1 = in1;
                k2 = in2;
            else
                ks = intersect(in1,in2);
                k1 = setdiff(in1,ks);
                k2 = setdiff(in2,ks);
            end
        end

        function [idxk1,idxks1,idxks2,idxk2] = getKeyIndices(in1,in2)
            [k1,ks,k2] = memZono.getUniqueKeys(in1,in2);

            % [~,idxk1] = ismember(k1,in1);
            % [~,idxk2] = ismember(k2,in2);
            idxk1 = zeros(1,length(k1));
            idxk2 = zeros(1,length(k2));
            for k = 1:length(k1)
                idxk1(k) = find(strcmp(in1,k1{k}));
            end
            for k = 1:length(k2)
                idxk2(k) = find(strcmp(in2,k2{k}));
            end
            if isempty(ks)
                idxks1 = []; idxks2 = [];
            else
                % [~,idxks1] = ismember(ks,in1);
                % [~,idxks2] = ismember(ks,in2);
                idxks1 = zeros(1,length(ks));
                idxks2 = zeros(1,length(ks));
                for k = 1:length(ks)
                    idxks1(k) = find(strcmp(in1,ks{k}));
                    idxks2(k) = find(strcmp(in2,ks{k}));
                end
            end
        end

    end
        

        %% Tags

        % function out = get.factorTags(obj); out = obj.tags.factors; end
        % function out = get.dimTags(obj); out = obj.tags.dims; end
        % function out = get.conTags(obj); out = obj.tags.cons; end


        % function obj = set.factorTags(obj,in)
        %     if isstruct(in)
        %         obj.tags.factors = in;
        %     else isalpha_num(in)
        %         if length(in) == obj.nG
        %             obj.tags.factors.keys = in;
        %         else
        %             obj.tags
        %     end
        % end



    methods
        %% Set Operations
        obj = minSum(obj1,obj2);
        obj = linMap(obj,M);
        obj = cartProd(obj1,obj2);
        obj = generalizedIntersection(obj1,obj2,R);
        obj = labeledIntersection(obj1,obj2,dims);



        %% Ploting
        plot(obj,dims,varargin);


        %% Overloading
        function obj = plus(in1,in2)
            % TODO: add if unlabeled
            obj = minSum(in1,in2);
        end
        function obj = mtimes(in1,in2)
            % TODO: assuming forward... include other as well
            if class(in2) == 'memZono'
                obj = in2.linMap(in1);
            else
                error('not coded')
            end
        end
        function obj = and(obj1,obj2)
            obj = labeledIntersection(obj1,obj2); %<-- intersects shared dims
        end
        % function obj = or(obj1,obj2)
        %     error('Union Not Coded')
        %     % obj = union(obj1,obj2);
        % end
                
        % Extended CartProd
        function obj = vertcat(varargin)
            obj = varargin{1};
            for i = 2:nargin
                obj = cartProd(obj,varargin{i});
            end
        end




        %% Indexing
        B = subsref(A,S);
        % A = subsasgn(A,S,B); %<---- not completed
        
        function out = projection(in,dims)
            if ~ismember(dims,in.dimKeys)
                dims = in.dimKeys(startsWith(in.dimKeys,dims));
            end
            [~,idx] = ismember(dims,in.dimKeys);
            keys_ = in.keys;
            keys_.dims = dims;
            out = memZono(in.G(idx,:),in.c(idx,:),in.A,in.b,in.vset,keys_);
        end



    end








end
