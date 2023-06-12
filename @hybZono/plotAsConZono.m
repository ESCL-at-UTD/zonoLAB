function [v,f] = plotAsConZono(obj,opts)

% Standardized header

[leaves] = getLeaves(obj,opts);
nLeaves = size(leaves,2);
opt = plotOptions('Display','off','SolverOpts',opts);
v = [];
f = [];
nVerts = zeros(nLeaves,1);
waitbarHandle = waitbar(0,['Plotting hybrid zonotope with ',num2str(nLeaves),' leaves.']);
for i = 1:nLeaves
    Zi = conZono(obj.Gc,obj.c+obj.Gb*leaves(:,i),obj.Ac,obj.b-obj.Ab*leaves(:,i));
    [vi,fi] = plot(Zi,opt);
    nVerts(i) = size(vi,1);
    v = [v;vi];
    if size(fi,2) > size(f,2)
        f = [f nan(size(f,1),size(fi,2)-size(f,2))]; 
    end
    if size(fi,2) < size(f,2)
        fi = [fi nan(1,size(f,2)-size(fi,2))];
    end
    f = [f;fi+sum(nVerts(1:i-1))];
    waitbar(i/nLeaves,waitbarHandle)
end
close(waitbarHandle)
% f = nan*ones(nLeaves,max(nVerts));
% count = 1;
% for i = 1:nLeaves
%     f(i,1:nVerts(i)) = count+[0:nVerts(i)-1];
%     count = count + nVerts(i);
% end

end