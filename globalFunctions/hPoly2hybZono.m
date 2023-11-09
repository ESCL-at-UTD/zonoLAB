% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Global function:
%       Returns a hybrid zonotope based on a collection of polytopes in H-rep
%   Syntax:
%       Z = hPoly2hybZono(H_collection)
%   Inputs:
%       H_collection - N x 1 cell array with i^th cell [H_i f_i] defining H-rep set H_i x <= f_i
%                      where H_i is a nH_i x n matrix and f_i is a nH_i x 1 vector
%   Outputs:
%       Z - hybrid zonotope in R^n
%   Notes:
%       Hybrid zonotope represents the union of the H-rep polytopes stored
%       in H_collection
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = hPoly2hybZono(H_collection)

N = length(H_collection);

out = hPoly2conZono(H_collection{1});
for i = 2:N
    out = union(out,hPoly2conZono(H_collection{i}));
end