clc;
clf;
clear;
n = 2; % 2 or 3
zono_type = 'hybZono';
% zono_type = 'conZono';

fprintf('Plot 1:\n')
% Zonotope with redundant constraints
subplot(2,2,1);
if strcmp(zono_type,'hybZono')
    Z = hybZono(eye(n),zeros(n,1),zeros(n,1),ones(2,n),ones(2,1),[0;0]);
else
    Z = conZono(eye(n),zeros(n,1),ones(2,n),[0;0]);
end
plot(Z,'b',0.8)
title('b = [0;0]')
axis(repmat([-1 1],1,n))

fprintf('------------------------\n')
fprintf('Plot 2:\n')
% Same zonotope but made infeasible by adjusting b above the tolerance level
% Perceived as infeasible
subplot(2,2,2);
Z.b(2) = 1e-4;
plot(Z,'b',0.8)
title('b = [0;1e-4]')
axis(repmat([-1 1],1,n))

fprintf('------------------------\n')
fprintf('Plot 3:\n')
% Same zonotope but made infeasible by adjusting b near the tolerance level
% One leaf is perceived as infeasible
subplot(2,2,3);
Z.b(2) = 1e-6;
[v,f]=plot(Z,'b',0.8);
title('b = [0;1e-6]')
axis(repmat([-1 1],1,n))

fprintf('------------------------\n')
fprintf('Plot 4:\n')
% Same zonotope but made infeasible by adjusting b below the tolerance level
% Perceived as redundant constraints
subplot(2,2,4);
Z.b(2) = 1e-7;
plot(Z,'b',0.8)
title('b = [0;1e-7]')
axis(repmat([-1 1],1,n))