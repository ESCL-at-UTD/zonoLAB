function B = subsref(A, S)
    switch S(1).type
        case {'.','{}'}
            B = builtin('subsref', A, S);
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
            i = getKeyIndices(i,A.dimKeys);
            j = getKeyIndices(j,A.factorKeys);
            k = getKeyIndices(k,A.conKeys);

            G_ = A.G(i,j);
            c_ = A.c(i,:);
            A_ = A.A(k,j);
            b_ = A.b(k,:);
            vset_ = A.vset(:,k);

            if ischar(i)
                keys_.dims = A.dimKeys;
            else
                keys_.dims = A.dimKeys(i);
            end
            if ischar(j)
                keys_.factors = A.factorKeys;
            else
                keys_.factors = A.factorKeys(j);
            end
            if ischar(k)
                keys_.cons = A.conKeys;
            else
                keys_.cons = A.conKeys(k);
            end

            B = memZono(G_,c_,A_,b_,vset_,keys_);
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