function [v,f] = plotHybZono1D(obj,opts)

% Standardized header


[v,f] = plotAsConZono(obj,opts);

% [leaves] = getLeaves(obj,opts);
% nLeaves = size(leaves,2);
% opt = plotOptions('Display','off');
% v = [];
% nVerts = zeros(nLeaves,1);
% for i = 1:nLeaves
%     Zi = conZono(obj.Gc,obj.c+obj.Gb*leaves(:,i),obj.Ac,obj.b-obj.Ab*leaves(:,i));
%     [vi,~] = plot(Zi,opt);
%     v = [v;vi];
%     nVerts(i) = size(vi,1);
% end
% f = nan*ones(nLeaves,max(nVerts));
% count = 1;
% for i = 1:nLeaves
%     f(i,1:nVerts(i)) = count+[0:nVerts(i)-1];
%     count = count + nVerts(i);
% end

end