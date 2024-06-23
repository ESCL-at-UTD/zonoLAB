% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Global function:
%       Returns a constrained zonotope based on a polytope in H-rep
%   Syntax:
%       Z = hPoly2conZono([H f])
%   Inputs:
%       H - nH x n matrix for set in H-rep (H x <= f)
%       f - nH x 1 vector for set in H-rep (H x <= f)
%   Outputs:
%       Z - constrained zonotope in R^n
%   Notes:
%       H-rep polytope of the form H x <= f.
%       Requires solving 2n linear programs.
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = hPoly2conZono(in)

H = in(:,1:end-1);
f = in(:,end);

n = size(H,2);

upperBounds = zeros(n,1);
lowerBounds = zeros(n,1);

for i = 1:n
    dir = zeros(1,n);
    dir(i) = 1;
    % Find upper bounds
    [x,~,~] = solveLP(dir,H,f,[],[],-realmax*ones(n,1),realmax*ones(n,1),[]);
    if isnan(x)
        error('H-rep is infeasible.')
    end
    upperBounds(i) = x(i);
    % Find lower bounds
    [x,~,~] = solveLP(-dir,H,f,[],[],-realmax*ones(n,1),realmax*ones(n,1),[]);
    if isnan(x)
        error('H-rep is infeasible.')
    end
    lowerBounds(i) = x(i);
end

G = diag(upperBounds-lowerBounds)/2;
c = (upperBounds+lowerBounds)/2;

out = zono(G,c);
out = halfspaceIntersection(out,H,f);

end