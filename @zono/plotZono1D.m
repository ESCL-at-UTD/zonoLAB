function [v,f] = plotZono1D(obj)

% Standardized header

v = obj.c + sum(abs(obj.G))*[-1 0; 1 0]; % Plot in x-direction with 0 in y-direction  
f = [1 2];

end