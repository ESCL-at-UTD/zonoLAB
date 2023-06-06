function plot(obj,varargin)

% Standardized header

opts = plotOptions;
if length(varargin) == 1
    if isa(varargin{1},'plotOptions')
        opts = varargin{1};
    else
        opts.FaceColor = varargin{1};
    end
elseif length(varargin) == 2
    opts.FaceColor = varargin{1};
    opts.FaceAlpha = varargin{2};
end
P = plotOptionsStruct(opts);

if (obj.n > 3) || (obj.n == 1)
    disp(['Can only plot in 2 or 3 dimensions.'])
elseif obj.n == 2
    obj.G(:,obj.G(2,:)<0) = -1*obj.G(:,obj.G(2,:)<0); % All generators in quadrants I and II
    angles = atan2(obj.G(2,:),obj.G(1,:));  % Compute angles
    [~,order] = sort(angles,'descend');     
    obj.G = obj.G(:,order);                 % Sort angles in clockwise order
    firstVertex = sum(obj.G');              % Compute first vertex
    % Compute remaining vertices by subtracting off one generator at a time
    vertsA = firstVertex - 2* [zeros(1,2);cumsum(obj.G,2)']; 
    vertsB = -vertsA; % Get centrally-symmetric vertices
    v = obj.c' + [vertsA;vertsB(2:end-1,:)]; % Shift by center
    P.XData = v(:,1);
    P.YData = v(:,2);
    patch(P)
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
    F = reshape([1:size(V,1)],4,[])';
    P.Faces = F;
    P.Vertices = obj.c' + V;
    patch(P)
end
    
end

