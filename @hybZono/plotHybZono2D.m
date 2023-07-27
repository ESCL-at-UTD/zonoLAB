% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return vertices and faces for a hybrid zonotope in 2D
%   Syntax:
%       [v,f] = plotHybZono2D(Z,optSolver)
%   Inputs:
%       Z - 2D hybrid zonotope in HCG-Rep (hybZono object)
%       optSolver - solver options needed for linear and mixed-integer linear propgrams
%   Outputs:
%       v - nV x 2 matrix, each row denoting the x (first column) and y (second column) positions
%                          of the nV vertices
%       f - nF x nMax matrix, each row denoting the vertices (up to nMax) contained
%                          in the nF faces (padded with NaN if face
%                          contains less than nMax vertices)
%   Notes:
%       Not intended to be called directly by user.
%       Use [v,f] = plot(obj,varargin) instead (method of abstractZono)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [v,f] = plotHybZono2D(obj,optSolver)

% Identify vertices and faces by treating hybrid zonotope at the union of
% constrained zonotopes
[v,f] = plotAsConZono(obj,optSolver);

end