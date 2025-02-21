function obj = orderReduction(obj,setType, orderReductionTechnique,final_order)
    % Default values for setType and orderReductionTechnique if not provided
if nargin < 3|| isempty(setType)
    setType = 'Zonotope';  % Default setType
end
if nargin < 4 || isempty(orderReductionTechnique)
    orderReductionTechnique = 'Exact';  % Default orderReductionTechnique
end

disp(['Set type: ', setType]);
disp(['Using order reduction technique: ', orderReductionTechnique]);

switch orderReductionTechnique
        case 'Exact'
            obj = exactApproximation(obj, setType);  % Call exact approximation
        case 'Outer'
            obj = outerApproximation(obj, setType,final_order);  % Call outer approximation
        case 'Inner'
            obj = innerApproximation(obj, setType,final_order);  % Call inner approximation
        otherwise
            error('Invalid reduction technique specified');
end
end

function obj= exactApproximation(obj,setType)
switch setType
    case 'Zonotope'
        nG = obj.nG;
        reduced_G = [];
        for i = 1:nG
            if any(obj.G(:,i))
                reduced_G = [reduced_G obj.G(:,i)];
            end
        end
        obj=zono(reduced_G,obj.c);
    case 'Constrained Zonotope'    
        A=obj.A;
        b=obj.b;
        c=obj.c;
        G=obj.G;
        while true
            [G,c,A,b,eff_r]=redundancy_remover_conzono(G,c,A,b);
            if isempty(eff_r)
                break
            end       
            i=length(eff_r);
            remove_position=eff_r(i);
            n_c=size(A,1);
            removable_constraint=0;
            for idx=1:n_c
                if A(idx,remove_position)~=0
                    removable_constraint=idx;
                    break
                end
            end
            if removable_constraint==0
            continue;
            end
            [G,c,A,b]=constraint_remover_function(G,c,A,b,remove_position,removable_constraint);
            obj=conZono(G,c,A,b);
        end
        obj=conZono(G,c,A,b);
        
    case 'Hybrid Zonotope'
        obj=obj
end
end

function obj= outerApproximation(obj,setType,final_order)
switch setType
    case 'Zonotope'
        G = obj.G;
        c=obj.c;  
        X = [G, -G]';  
        Co = X' * X;   
        % Step 4: Perform Singular Value Decomposition (SVD)
        [U, ~, ~] = svd(Co);  % U contains the principal components (eigenvectors)
     
        Z_transformed=U'*obj;    
        obj=U*boundingBox(Z_transformed);
    case 'Constrained Zonotope'
        obj=obj
    case 'Hybrid Zonotope'
        obj=obj
end
end

function obj= innerApproximation(obj,setType,final_order)
switch setType
    case 'Zonotope'
        G=obj.G;
        c=obj.c;
        n_g=size(G,2);
        n_r=final_order;
        
        %% Sorting generators and computing T
        % Compute generator norms
        Norms = vecnorm(G);
        % Sort generator norms in decreasing order
        [Norms_sort,indx_sort] = sort(Norms,'descend');
        G = G(:,indx_sort);
        
        % Seperate n_r largest generators
        Z1.G = G(:,1:n_r);
        Z2.G = G(:,n_r+1:end);
        
        % Compute magnitude of dot product between generators in Z1 and Z2
        alpha_abs = zeros(n_r,n_g-n_r);
        alpha = zeros(n_r,n_g-n_r);
        for i = 1:n_r
            for j = 1:n_g-n_r
                alpha_abs(i,j) = abs(Z1.G(:,i)'*Z2.G(:,j));
                alpha(i,j) = Z1.G(:,i)'*Z2.G(:,j);
            end
        end
        % Normalize the dot product with respect to largest dot product
        alpha_norm = alpha*diag(1./max(alpha_abs,[],1));
        
        T2 = zeros(n_r,n_g-n_r);
        T2(alpha_norm==1) = 1;
        T2(alpha_norm==-1) = -1;
        
        T = [eye(n_r);T2'];
        G_transformed=G*T;
        obj=zono(G_transformed,c);
end
end