% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Helper Function for Zono2HPoly that removes columns from M
%       specified in List
%   Example: M = [1 2 3 4 5;          List = [2 5]
%                 6 7 8 9 0];
%           out = [1 3 4;
%                  6 8 9];
%   Syntax:
%       out = RemoveColumns(M, List)
%   Inputs:
%       M - n by m matrix
%       List - 1 by p matrix. Assumed to be sorted in ascending order (e.g.
%           List = [1 2 7 9], not List = [2 9 1 7] etc.).
%   Outputs:
%       out - n by (m-p) matrix
%
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function out = removeColumns(M, List)
% Need input error checking here
%  1) Check that List is a row vector
%  2) Run a sort on List
%  3) Check that all entires in List are within size(M,2)

    p = size(List,2);
    out = M;
    for k = 1:p
        out(:, List(k) - (k-1)) = [];
    end
end