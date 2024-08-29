% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%   Method:
%       Returns the convex hull of two zonotopic sets, Z = CH(X \cup Y), or
%           the convex hull of a single zonotopic set. 
%   Syntax:
%       Z = convexHull(X,Y)
%   Inputs:
%       X - zonotopic set in R^n (hybZono, conZono, or zono object)
%       Y - zonotopic set in R^n (hybZono, conZono, or zono object)
%   Outputs:
%       Z - zonotopic set in R^n (conZono object)
%   Notes:
%       As currently written, can only perform convex hull of HZ in 1, 2,
%       or 3 dimensions (as it relies on the plotting function)
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
function out = convexHull(obj1,varargin)

switch nargin
    case 2
        obj2 = varargin{1};
    case 1
        if isa(obj1, 'zono')
            out = conZono(obj1);
        elseif isa(obj1, 'conZono')
            out = obj1;
        elseif isa(obj1, 'hybZono')
            out = INTERNAL_fcn_hybZono_hull(obj1);
        end
        return
    otherwise
        error('The method convexHull accepts only one or two inputs.')
end

        if isa(obj1,'zono')
            obj1 = conZono(obj1);
        end
        if isa(obj2,'zono')
            obj2 = conZono(obj2);
        end
        if isa(obj1, 'hybZono')
            obj1 = INTERNAL_fcn_hybZono_hull(obj1);
        end
        if isa(obj2, 'hybZono')
            obj2 = INTERNAL_fcn_hybZono_hull(obj2);
        end

        c = (obj1.c + obj2.c)/2;
        G = [obj1.G obj2.G (obj1.c-obj2.c)/2 zeros(obj1.n,2*(obj1.nG+obj2.nG))];

        H = [eye(obj1.nG) zeros(obj1.nG,obj2.nG)  -0.5*ones(obj1.nG,1);...
            -eye(obj1.nG) zeros(obj1.nG,obj2.nG)  -0.5*ones(obj1.nG,1);...
             zeros(obj2.nG,obj1.nG)  eye(obj2.nG)  0.5*ones(obj2.nG,1);...
             zeros(obj2.nG,obj1.nG) -eye(obj2.nG)  0.5*ones(obj2.nG,1)];
        I = eye(2*(obj1.nG + obj2.nG));
        f = -0.5*ones(2*(obj1.nG + obj2.nG),1);

        A = [obj1.A zeros(obj1.nC,obj2.nG) -obj1.b/2 zeros(obj1.nC,size(G,2)-obj1.nG-obj2.nG-1);...
            zeros(obj2.nC,obj1.nG) obj2.A   obj2.b/2 zeros(obj2.nC,size(G,2)-obj2.nG-obj1.nG-1);...
            H I];
        b = [obj1.b/2;obj2.b/2;f];

        out = conZono(G,c,A,b);
end

function out = INTERNAL_fcn_hybZono_hull(Z)
    [v,~] = plot(Z, plotOptions('Display','off'));
    v_conv = v(convhull(v),:);
    out = vPoly2Zono(v_conv');

    % leaves = getLeaves(Z, solverOptions);
    % if isempty(leaves)
    %     warning('zonoLAB:EmptyZonotope','Hybrid zonotope is empty and cannot be plotted.')
    %     out = conZono([],[],[],[]);
    %     return
    % end
    % nLeaves = size(leaves,2);
    % out = conZono(Z.Gc, Z.c + Z.Gb*leaves(:,1), Z.Ac, Z.b - Z.Ab*leaves(:,1));
    % ax = gca;
    % plot(out, 'm',1);
    % if nLeaves > 1
    %     for i = 2:nLeaves
    %         Zi = conZono(Z.Gc, Z.c + Z.Gb*leaves(:,i), Z.Ac, Z.b - Z.Ab*leaves(:,i));
    %         if i>2
    %             ax.Children(2).delete
    %         end
    %         plot(Zi, 'g', 1);
    %         out = convexHull(out, Zi)
    %         ax.Children(2).delete
    %         plot(out, 'm', 0.4);
    %     end
    % end

end