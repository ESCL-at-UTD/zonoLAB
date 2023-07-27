% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return vertices and faces for a zonotope in 1D
%   Syntax:
%       [v,f] = plotZono1D(Z)
%   Inputs:
%       Z - 1D zonotope in G-Rep (zono object)
%   Outputs:
%       v - 2 x 2 matrix, first row indicating lower and upper bounds (x-direction),
%                         second row indicting zeros for plotting purposes (y-direction)
%       f - 1 x 2 vector, indicating a single face containing both vertices 
%   Notes:
%       Not intended to be called directly by user.
%       Use [v,f] = plot(obj,varargin) instead (method of abstractZono)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [v,f] = plotZono1D(obj)

% Assume zonotope spans the x-direction (horizontal) with 0 in the y-direction
v = obj.c + sum(abs(obj.G))*[-1 0; 1 0];    % Upper and lower bounds
f = [1 2];                                  % Face

end