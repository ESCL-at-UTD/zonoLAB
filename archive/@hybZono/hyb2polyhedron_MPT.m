%	Hybrid zonotope object function: convert a hybrid zonotope to an
%		equivalent complex of convex polyhedron using MPT 3.0
% 
%	hybZono/hyb2polyhedron_MPT
% 
%	Syntax: 
%		P = hyb2polyhedron_MPT(Z,combos)
% 
%	Inputs:
%		Z : n dimensional hybrid zonotope object in HCG-rep
% 		combos: N feasible combinations of the binary factors (optional, if
%			not provided function will call getLeaves to find them)
% 
%	Outputs:
% 		P : complex of N convex polyhedron, one for each feasible
%			combination of binary factors
%	
% Notes: 
% 
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function P = hyb2polyhedron_MPT(obj,combos)
Z = obj;
Z = getDimensions(Z);

% find possible combination of binaries then convert to convex of polyhedron using MPT
if Z.ngb
	
	if nargin < 2
		combos = getLeaves(Z);
	end
	
	% convert the hybrid zonotpe to an array of polyhedron using the possible combinations of binary factors
	for i = 1:size(combos,2)
		% find continuous constrained hybercube to be projected
		if Z.nc
			Box = Polyhedron('lb',-ones(Z.ngc,1),'ub',ones(Z.ngc,1),'He',[Z.Ac Z.b-Z.Ab*combos(:,i)]);
		else
			Box = Polyhedron('lb',-ones(Z.ngc,1),'ub',ones(Z.ngc,1));
		end
		if Box.isEmptySet
			warning('box is empty')
		end
		% find affine image of constrained hypercube with binary factors
		mappedBox = affineMap(Box,Z.Gc);
		Zj = plus(Z.c+Z.Gb*combos(:,i),mappedBox);
		if Zj.isEmptySet
			warning('mapped Polyhedron is empty')
		end
		if i > 1
			P = [ P , Zj ];
		else
			P = Zj;
		end
	end

else

	% this is just a constrained zonotope
	% find continuous constrained hybercube to be projected
	if Z.nc
		Box = Polyhedron('lb',-ones(Z.ngc,1),'ub',ones(Z.ngc,1),'He',[Z.Ac Z.b]);
	else
		Box = Polyhedron('lb',-ones(Z.ngc,1),'ub',ones(Z.ngc,1));
	end
	% find affine image of constrained hypercube
	P = plus(Z.c,affineMap(Box,Z.Gc));

end
	
end % end hyb2polyhedron