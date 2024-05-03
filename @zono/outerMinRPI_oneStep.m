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
    nw = W.nG;
    ng = nw*(s+1);

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

    % Setup Optimization Problem
    prob = optimproblem;
    prob.Constraints.gamma_1 = A*G*Phi == G*Gamma_1;
    prob.Constraints.gamma_2 = W.G == G*Gamma_2;
    prob.Constraints.beta = (eye(n) - A)*W.c - c == G*beta;

    
    % abs() linear constraints
    g = optimvar('g',ng,(ng+nw));
    b = optimvar('b',ng);
    Gamma = horzcat(Gamma_1,Gamma_2);
    prob.Constraints.g = vertcat(g,g) >= vertcat(Gamma,-Gamma);
    prob.Constraints.b = vertcat(b,b) >= vertcat(beta,-beta);
    prob.Constraints.sum = horzcat(g,b)*ones(ng+nw+1,1)<= Phi*ones(ng,1);

    
    % objective version
    phi_max = optimvar('phi_max',1);
    prob.Constraints.phi_max = phi <= phi_max*ones(ng,1);
    prob.Objective = phi_max;%<== inf norm since $phi_{i} > 0 \forall_{i}$

    % Solve for result
    sol = solve(prob);

    F_approx = zono(G*diag(sol.phi),sol.c);

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