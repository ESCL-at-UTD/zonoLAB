function varargout = subsref(obj, S)
    % Overload of varargout: https://www.mathworks.com/help/matlab/matlab_oop/code-patterns-for-subsref-and-subsasgn-methods.html
    switch S(1).type
        case {'.','{}'}
            [varargout{1:nargout}] = builtin('subsref', obj, S);
        case '()'
            switch numel(S(1).subs)
                case 1
                    i = S(1).subs{1}; j = ':'; k = ':';
                case 2
                    i = S(1).subs{1}; j = S(1).subs{2}; k = ':';
                case 3
                    i = S(1).subs{1}; j = S(1).subs{2}; k = S(1).subs{3};
                otherwise
                    error('indexing not specificed')
            end

            % TODO: Need error handling here if key doesn't exist
            i = getKeyIndices(i,obj.dimKeys);
            j = getKeyIndices(j,obj.factorKeys);
            k = getKeyIndices(k,obj.conKeys);

            G_ = obj.G(i,j);
            c_ = obj.c(i,:);
            A_ = obj.A(k,j);
            b_ = obj.b(k,:);
            vset_ = obj.vset(:,k);

            if ischar(i)
                keys_.dims = obj.dimKeys;
            else
                keys_.dims = obj.dimKeys(i);
            end
            if ischar(j)
                keys_.factors = obj.factorKeys;
            else
                keys_.factors = obj.factorKeys(j);
            end
            if ischar(k)
                keys_.cons = obj.conKeys;
            else
                keys_.cons = obj.conKeys(k);
            end

            varargout{1} = memZono(G_,c_,A_,b_,vset_,keys_);
    end
end


function idx = getKeyIndices(in,keys)
    if isnumeric(in)
        idx = in;
    else
        if iscell(in)
            [~,idx] = ismember(in,keys);
        else
            if ~strcmp(in,':')
                [~,idx] = ismember(in,keys);
            else
                idx = ':';
            end
        end
    end
end