%	Hybrid zonotope object function: get / update the dimension properties 
%									of a hybrid zonotope
%		hybZono/getDimensions
% 
%	Syntax: 
%		Zh = getDimensions(Zh)
% 
%	Inputs:
%		Zh : n dimensional hybrid zonotope object in HCG-rep
% 
%	Outputs:
% 		Zh : n dimensional hybrid zonotope object in HCG-rep with updated 
%				dimension properties
% 
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function obj = getDimensions(obj)

obj.n  = size(obj.c,1);

if isempty(obj.Gc)
    obj.ngc = 0;
    obj.Gc = zeros(obj.n,obj.ngc);
else
    obj.ngc = size(obj.Gc,2);
end
if isempty(obj.Gb)
    obj.ngb = 0;
    obj.Gb = zeros(obj.n,obj.ngb);
else
    obj.ngb = size(obj.Gb,2);
end    
if isempty(obj.b)
    obj.nc = 0;
    obj.b = zeros(0,1);
else
    obj.nc = size(obj.b,1);
end
if isempty(obj.Ac)
    obj.Ac = zeros(obj.nc,obj.ngc);
end
if isempty(obj.Ab)
    obj.Ab = zeros(obj.nc,obj.ngb);
end

end