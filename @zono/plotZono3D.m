% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Return vertices and faces for a zonotope in 3D
%   Syntax:
%       [v,f] = plotZono3D(Z)
%   Inputs:
%       Z - 3D zonotope in G-Rep (zono object)
%   Outputs:
%       v - nV x 3 matrix, each row denoting the x (first column), y (second column), 
%                          and z (third column) positions of the nV vertices
%       f - nF x 4 matrix, each row denoting the four vertices contained
%                          in the nF faces 
%   Notes:
%       Not intended to be called directly by user.
%       Use [v,f] = plot(obj,varargin) instead (method of abstractZono)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function [v,f] = plotZono3D(obj)

if rank(obj.G) == 1
    nullVec = null(obj.G');
    crossVec = cross(nullVec(:,1),nullVec(:,2));
    reducedG = crossVec'*obj.G;
    reducedObj = zono(reducedG,0);
    optPlot = plotOptions('Display','off');
    [reducedV,~] = plot(reducedObj,optPlot);
    v = obj.c' + reducedV(:,1).*crossVec';
    f = [1 2];
elseif rank(obj.G) == 2
    reducedG = orth(obj.G)'*obj.G;
    reducedObj = zono(reducedG,zeros(2,1));
    optPlot = plotOptions('Display','off');
    [reducedV,~] = plot(reducedObj,optPlot);
    v = obj.c' + reducedV*orth(obj.G)';
    f = [1:size(v,1)];
else
    twoCombos = ff2n(2); % Identify all of the combinations of {-1,1}^2
    twoCombos(twoCombos == 0) = -1; 
    twoCombos(4:-1:3,:) = twoCombos(3:4,:); % Reorder so that facets plot correctly
    allCombos = nchoosek(1:obj.nG,2); % Find all combinations of two generators
    V = zeros(8*size(allCombos,1),3); % Initialize matrix to store vertices
    for i = 1:size(allCombos,1)
        verts = twoCombos*obj.G(:,allCombos(i,:))'; % Compute the four vertices for this facet
        FN = cross(obj.G(:,allCombos(i,1)),obj.G(:,allCombos(i,2))); % Find the vector normal to this facet
        remainingGens = setdiff([1:obj.nG],allCombos(i,:)); % Find the indices of the other generators 
        dotProds = FN'*obj.G(:,remainingGens); % Compute the dot product of these other generators with the facet normal vector
        signs = sign(dotProds); % Determine the sign of the resulting dot product
        centerA = signs*obj.G(:,remainingGens)'; % Position the center of the facet based on the other generators
        centerB = -signs*obj.G(:,remainingGens)'; % Centrally-symmetric
        V([1:4]+(i-1)*8,:) = centerA + verts; % Add vertices to list
        V([5:8]+(i-1)*8,:) = centerB + verts;
    end
    f = reshape([1:size(V,1)],4,[])';
    v = obj.c' + V;
end

end