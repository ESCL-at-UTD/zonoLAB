function A = subsasgn(A, S, B)
    switch S(1).type
        case {'.','{}'}
            A = builtin('subsasgn', A, S, B);
        case '()'
            switch numel(S(1).subs)
                case 1
                    i = S(1).subs{1}; j = ':'; k = ':';
                case 2
                    i = S(1).subs{1}; j = S(1).subs{2}; k = ':';
                case 3
                    i = S(1).subs{1}; j = S(1).subs{2}; k = ':';
                otherwise
                    error('indexing not specificed')
            end

            %% conZonoM version 
            % % ensure numeric
            % if isnumeric(i); i = A.dimKeys(i); end
            % if isnumeric(j); j = A.factorKeys(j); end
            % if isnumeric(k); k = A.conKeys(k); end
            % 
            % % Asign to values
            % A.c_dict(i,1) = B.c_dict;
            % A.G_dict(i,j) = B.G_dict;
            % A.A_dict(k,j) = B.A_dict;
            % A.b_dict(k,1) = B.b_dict;

            warning("this asignment is untested... probably won't work")
    end
end