% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Helper Function for Zono2HPoly that computes the n-dimensional
%       cross-product of a matrix H
%   Syntax:
%       out = nCross(H)
%   Inputs:
%       H - n by n-1 matrix
%   Outputs:
%       out - n by 1 vector
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function out = nCross(H)
    if length(size(H)) ~= 2
        error('Input matrix H must be two-dimensional.')
    end
    if size(H,1) ~= size(H,2)+1
        error('Input matrix H must have dimensions n by n-1.')
    end

    n = size(H,1);
    out = NaN(n,1);
    for i = 1:n
        H_i = [H(1:i-1, :); H(i+1:n, :)];
        out(i) = ((-1)^(i+1))*det(H_i);
    end
end