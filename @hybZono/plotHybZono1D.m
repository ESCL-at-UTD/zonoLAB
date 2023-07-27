% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return vertices and faces for a hybrid zonotope in 1D
%   Syntax:
%       [v,f] = plotHybZono1D(Z,optSolver)
%   Inputs:
%       Z - 1D hybrid zonotope in HCG-Rep (hybZono object)
%       optSolver - solver options needed for linear and mixed-integer linear propgrams
%   Outputs:
%       v - nV x 2 matrix, each row denoting the x (first column) and y (second column) positions
%                          of the nV vertices
%       f - nF x 2 matrix, each row denoting the two vertices contained
%                          in the nF faces
%   Notes:
%       Not intended to be called directly by user.
%       Use [v,f] = plot(obj,varargin) instead (method of abstractZono)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [v,f] = plotHybZono1D(obj,optSolver)

% Identify vertices and faces by treating hybrid zonotope at the union of
% constrained zonotopes
[v,f] = plotAsConZono(obj,optSolver);

end