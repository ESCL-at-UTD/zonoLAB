%	Hybrid zonotope object function: remove zero rows and columns
% 
%	hybZono/removeZeros
% 
%	Syntax: 
%		Zh = removeZeros(Zh,tol)
% 
%	Inputs:
%		Zh : n dimensional hybrid zonotope in HCG-rep
%		tol : (optional) scalar with tolerance of what is considered to be 
%				zero default is tol = 1e-12
% 
%	Outputs:
% 		Zh : n dimensional hybrid zonotope in HCG-rep without zeros rows or columns
%	
%	Trevor Bird - bird6@purdue.edu - Purdue 2022
%   Modified: Jonah Glunt, jjg57@psu.edu
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function Zh = removeZeros(Zh,tol)

if nargin < 2
	tol = 1e-12;
end

% are there zeros to remove?
if isempty(Zh.Ac) && isempty(Zh.Ab) && isempty(Zh.b)
	rowZero = [];
else
	rowZero = find(all(abs([ Zh.Ac Zh.Ab Zh.b ])<tol,2));
end
colZeroC = find(all(abs([ Zh.Gc ; Zh.Ac  ])<tol,1));
colZeroB = find(all(abs([ Zh.Gb ; Zh.Ab  ])<tol,1));

% remove empty constraints
Zh.Ac(rowZero,:) = [];
if ~isempty(Zh.Ab)
	Zh.Ab(rowZero,:) = [];
end
Zh.b(rowZero) = [];

% remove zero columns
Zh.Gc(:,colZeroC) = [];
if ~isempty(Zh.Ac)
	Zh.Ac(:,colZeroC) = [];
end
Zh.Gb(:,colZeroB) = [];
if ~isempty(Zh.Ab)
	Zh.Ab(:,colZeroB) = [];
end

% % if we don't have binary generators anymore reset binTree
% if ~Zh.nGb
% 	Zh.binTree = [];
% end

end