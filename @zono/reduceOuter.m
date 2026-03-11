function obj = reduceOuter(obj)
% REDUCEOUTER Computes an outer approximation of a zonotope
% using principal directions and a bounding box.

    G = obj.G;
    c = obj.c;

    % Construct symmetric generator set
    X = [G -G]';

    % Covariance like matrix
    Co = X' * X;

    % Principal directions
    [U,~,~] = svd(Co);

    % Rotate zonotope into principal basis
    Z_transformed = U' * obj;

    % Compute bounding box in rotated space
    Z_box = boundingBox(Z_transformed);

    % Transform back
    obj = U * Z_box;

end