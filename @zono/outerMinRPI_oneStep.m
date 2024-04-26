% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Computes an outer approximation of the minimal Robust Positively 
%       Invariant (mRPI) set for an autonomous, discrete-time system under 
%       a bounded disturbance set.
%       Based on Theorem 6 from Vignesh2022
%   Syntax:
%       F_approx = outerMinRPI_oneStep(A, W, s)
%   Inputs:
%       A - square, Schur system matrix
%       W - zonotope representation of the bounded disturbance set,
%           centered at the origin. 
%           Function does not currently support conZono or hybZono W.
%       s - Describes the highest power of A used in the generator matrix
%   Outputs:
%       F_approx - outer-approximated mRPI set as a zonotope
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [F_approx] = outerMinRPI_oneStep(A, W, s, options)
    arguments
        A {mustBeSquareMatrix,mustBeSchurSystemMatrix}
        W zono
        s (1,1) {mustBeInteger} = 5;
    end
    % options
    arguments
        options.verbose = false;
    end


    
    % sizes
    n = size(A,1);
    if W.n ~= n; error('Input set should be same size as update matrix'); end
    ng = n*(s+1);
    nw = W.nG;

    % Template Generator Set
    for i = 1:s+1
        G_{i} = A^(i-1) * W.G;
    end
    G = horzcat(G_{:});


    % Optimvar setup
    phi = optimvar('phi',ng,'LowerBound',0);
    Phi = diag(phi);
    c = optimvar('c',n);
    Gamma_1 = optimvar('gamma_1',ng,ng);
    Gamma_2 = optimvar('gamma_2',ng,nw);
    beta = optimvar('beta',ng);

    prob = optimproblem;
    prob.Objective = max(phi);%<== inf norm since $phi_{i} > 0 \forall_{i}$

    prob.Constraints.gamma_1 = A*G*Phi == G*Gamma_1;
    prob.Constraints.gamma_2 = W.G == G*Gamma_2;
    prob.Constraints.beta = (eye(n) - A)*W.c - c == G*beta;

    Gamma = horzcat(Gamma_1,Gamma_2);

    % abs() linear constraints
    g = optimvar('g',ng,(ng+nw));
    b = optimvar('b',n_g);

    prob.Constraints.g = vertcat(g,g) >= vertcat(Gamma,-Gamma);
    prob.Constraints.b = vertcat(b,b) >= vertcat(beta,-beta);
    prob.Constraints.sum = g*ones(ng+nw,1) + b <= Phi*ones(ng);

    sol = solve(prob);

    F_approx = zono(G*diag(sol.phi),sol.c);




    
    % Input Error Checking

    %d = abs(eig(A));
    %for i = 1:length(d)
    %    if d > 1
    %        error('Input A must be a discrete-time stable (Schur) matrix.')
    %    end
    %end

    % %Check if W is a zonotope
    % if class(W) ~= 'zono'
    %     error('Input set W must be a zonotope object.');
    % end

    % %Check if epsilon is a scalar
    % if (~isscalar(epsilon))
    %     error('Input epsilon must be a scalar.')
    % end

    % verbose = false;
    % if ~isempty(varargin)
    %     verbose = true;
    % end

    % %RPI Algorithm Start
    % s = 0;
    % [H, f] = zono2HPoly(W);
    % nH = size(H,1);
    % Fs = W;
    
    % for i = 1:max_iter
    %     s = s + 1;
    %     AsW = (A^s)*W;
    %     alpha_o = NaN(nH, 1);
        
    %     % Find alpha
    %     for j = 1:nH
    %         alpha_o(j) = supportFunc(AsW, H(j,:)')/f(j);
    %     end
    %     alpha = max(alpha_o);
        
    %     % Find M(s)
    %     Box = boundingBox(Fs);
    %     M_s = max(max(Box.G));
        
    %     % Check if alpha <= epsilon / (epsilon + M(s))
    %     if alpha <= (epsilon / (epsilon + M_s))
    %         break;
    %     else
    %         Fs = Fs + ((A^s)*W);
    %     end
    % end
    % if i >= max_iter
    %     disp('outerMinRPI could not converge to a minRPI approximation in the specified number of iterations.');
    % end
    
    % if i < max_iter && verbose
    %     str = 'minRPI approximation had successfully converged. s = %d, alpha = %f\n';
    %     fprintf(str,s,alpha);
    % end
    
    % F_approx = (1/(1-alpha))*Fs;
end



function mustBeSquareMatrix(A)
    sz = size(A);
    if (length(sz) > 2 || length(sz) < 2)
        error('Input A must be a two-dimensional matrix.')
    elseif sz(1) ~= sz(2)
        error('Input A must be a square matrix.')
    end
end

function mustBeSchurSystemMatrix(A)
    if any(abs(eig(A)) > 1)
        error('All eigen values must be stable') 
    elseif sum(abs(eig(A))==1) > 1 %<== does not check if repeated roots are the same place...
        error("Can't have repeated roots on the boundary")
    end
end