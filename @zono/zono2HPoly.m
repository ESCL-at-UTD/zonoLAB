% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the equivalent halfspace representation for a given 
%       zonotope Z 
%   Syntax:
%       [C, d] = Zono2HPoly(Z)
%   Inputs:
%       Z - zonotope object with generator matrix G (n by ng dimensions)
%       and center c (n by dimensions)
%   Outputs:
%       C = Matrix containing the normals of each halfspace
%       d = Vector containing the offsets of each halfspace
%   C and d are such that C*x <= d is an equivalent set to Z   
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [C, d] = zono2HPoly(obj, varargin)
    if class(obj) ~= 'zono'
        error('Function zono2HPoly only supports input objects of class zono.');
    end

    G = obj.G;
    c = obj.c;
    n = obj.n;
    ng = obj.nG;
    
    Remove_Column_List = nchoosek(1:ng, ng-(n-1));
    nH = nchoosek(ng, n-1);
    
    C_plus = NaN(nH, n);
    d_plus = NaN(nH, 1);
    d_minus = NaN(nH, 1);
    
    parallel = false;
    if ~isempty(varargin)
        parallel = true;
    end
    
    if parallel
        if nH <= 1000000
            disp('Parallel computation in zono2HPoly is only recommended for high values of nG Choose (n-1).');
            disp('Avoiding parallel computation may be faster for your computation.');
        end
        parfor i = 1:nH
            G_removed = removeColumns(G,Remove_Column_List(i,:));
            nX = nCross(G_removed);
            C_plus(i,:) = nX' / norm(nX,2);
            Delta_d_i = sum( abs(C_plus(i,:)*G ) );
            d_plus(i) = C_plus(i,:)*c + Delta_d_i;
            d_minus(i) = -C_plus(i,:)*c + Delta_d_i;
        end
    else
        if nH > 1000000
            disp('zono2HPoly has been called for a high value of nG Choose (n-1).');
            disp('Consider using the parallel option for faster computation time.');
        end
        for i = 1:nH
            G_removed = removeColumns(G,Remove_Column_List(i,:));
            nX = nCross(G_removed);
            C_plus(i,:) = nX' / norm(nX,2);
            Delta_d_i = sum( abs(C_plus(i,:)*G ) );
            d_plus(i) = C_plus(i,:)*c + Delta_d_i;
            d_minus(i) = -C_plus(i,:)*c + Delta_d_i;
        end
    end
    
    C = [C_plus; -C_plus];
    d = [d_plus; d_minus];
end