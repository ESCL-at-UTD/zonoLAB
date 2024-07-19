% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Computes an outer approximation of the minimal Robust Positively 
%       Invariant (mRPI) set for an autonomous, discrete-time system under 
%       a bounded disturbance set.
%       Based on Algorithm 1 from "Invariant Approximations of the Minimal
%       Robust Positively Invariant Set" (Rakovic et. al., 2005)
%   Syntax:
%       [F_approx, s, alpha] = Approx_RPI_Zono(A, W, epsilon, max_iter)
%   Inputs:
%       A - square, Schur system matrix
%       W - zonotope representation of the bounded disturbance set,
%           centered at the origin. 
%           Function does not currently support conZono or hybZono W.
%       epsilon - scalar, specifies the level of approximation accuracy
%       max_iter - scalar, specifies the max # of computation iterations. 
%       verbose (optional) - if true prints a statement about the result; default is false
%   Outputs:
%       F_approx - epsilon-close outer-approximated mRPI set as a zonotope
%       s - Describes the highest power of A used in the finite sum of F_approx (s-1)
%       alpha - Describes the scaling factor of 1/(1 - alpha) used on the
%       finite zonotope sum to reach F_approx.
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [F_approx, s, alpha] = outerMinRPI(A, W, epsilon, max_iter, varargin)
    % Input Error Checking
    sz = size(A);
    if (length(sz) > 2 || length(sz) < 2)
        error('Input A must be a two-dimensional matrix.')
    elseif sz(1) ~= sz(2)
        error('Input A must be a square matrix.')
    end

    %d = abs(eig(A));
    %for i = 1:length(d)
    %    if d > 1
    %        error('Input A must be a discrete-time stable (Schur) matrix.')
    %    end
    %end

    %Check if W is a zonotope
    if class(W) ~= 'zono'
        error('Input set W must be a zonotope object.');
    end

    %Check if epsilon is a scalar
    if (~isscalar(epsilon))
        error('Input epsilon must be a scalar.')
    end

    verbose = false;
    if ~isempty(varargin)
        verbose = true;
    end

    %RPI Algorithm Start
    s = 0;
    [H, f] = zono2HPoly(W);
    nH = size(H,1);
    Fs = W;
    
    for i = 1:max_iter
        s = s + 1;
        AsW = (A^s)*W;
        alpha_o = NaN(nH, 1);
        
        % Find alpha
        for j = 1:nH
            alpha_o(j) = supportFunc(AsW, H(j,:)')/f(j);
        end
        alpha = max(alpha_o);
        
        % Find M(s)
        Box = boundingBox(Fs);
        M_s = max(max(Box.G));
        
        % Check if alpha <= epsilon / (epsilon + M(s))
        if alpha <= (epsilon / (epsilon + M_s))
            break;
        else
            Fs = Fs + ((A^s)*W);
        end
    end
    if i >= max_iter
        disp('outerMinRPI could not converge to a minRPI approximation in the specified number of iterations.');
    end
    
    if i < max_iter && verbose
        str = 'minRPI approximation had successfully converged. s = %d, alpha = %f\n';
        fprintf(str,s,alpha);
    end
    
    F_approx = (1/(1-alpha))*Fs;
end