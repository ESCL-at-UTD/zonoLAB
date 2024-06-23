% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   This function will add the zonoLab toolbox to your MATLAB path.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function zonoLab_Install

addpath(cd)
addpath('examples');
addpath(genpath('examples'))
addpath('globalFunctions');
addpath(genpath('globalFunctions'))
addpath('settings');
addpath(genpath('settings'))

if exist('gurobi','dir') == 0
    disp('WARNING!')
    disp('It looks like you do not have GUROBI installed for MATLAB.')
    disp('zonoLab functionality is severely limited without GUROBI.')
end

end

