function plot(obj,varargin)
% Plot a constrained zonotope
% Default color: red. Default alpha: 1.
if isempty(varargin)
    color = 'r'; alpha = 1;
elseif length(varargin) == 2
    color = varargin{1}; alpha = varargin{2};
end

if obj.n > 3
    disp(['Cannot plot sets in ',num2str(obj.n),' dimensions!'])
else
    Box = Polyhedron('H',[obj.A obj.b]);
    P = plus(obj.c,affineMap(Box,obj.G));
    plot(P,'color',color,'alpha',alpha)
end
end