function [v,f] = plotZono2D(obj)

% Standardized header

if rank(obj.G) == 1
    nullVec = null(obj.G');
    reducedG = [-nullVec(2) nullVec(1)]*obj.G;
    reducedObj = zono(reducedG,0);
    opt = plotOptions('Display','off');
    [reducedV,~] = plot(reducedObj,opt);
    v = obj.c' + reducedV(:,1).*[nullVec(1) -nullVec(2)];
    f = [1 2];
else
    obj.G(:,obj.G(2,:)<0) = -1*obj.G(:,obj.G(2,:)<0); % All generators in quadrants I and II
    angles = atan2(obj.G(2,:),obj.G(1,:));  % Compute angles
    [~,order] = sort(angles,'descend');
    obj.G = obj.G(:,order);                 % Sort angles in clockwise order
    firstVertex = sum(obj.G');              % Compute first vertex
    % Compute remaining vertices by subtracting off one generator at a time
    vertsA = firstVertex - 2* [zeros(1,2);cumsum(obj.G,2)'];
    vertsB = -vertsA; % Get centrally-symmetric vertices
    v = obj.c' + [vertsA;vertsB(2:end-1,:)]; % Shift by center
    f = [1:size(v,1)];
end

end