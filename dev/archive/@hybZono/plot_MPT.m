%	Hybrid zonotope object function: plot a hybrid zonotope
%	
%	hybZono/plot
% 
%	Syntax: 
%		[ p , times ] = plot(Z,color,alpha,res,dims)
% 
%	Inputs:
%		Z : n dimensional hybrid zonotope object in HCG-rep
% 		color : [ r b g ] color to plot
%		alpha : scalar for the shadding of the plot
%		res : resolution for plotting inner approximation by solving res 
%				number of LPs to find max / min around unit half circle. 
%				DEFAULT : res == 0 => use MPT to convert from HCG to H -rep
%		dims : [ z1 , z2 ] dimensions of the set to be projected (optional
%			   only required if set is greater than 3D)
% 
%	Outputs:
% 		p : handle for plot 
%		times : struct with how long it took to do the three parts of
%			plotting -> .combos .convert .plot .all
%	
% Notes: If Z does not have the binTree already found this function will
% call getLeaves to find the nonempty leaves of the binary tree 
% 
% ToDo:
% Add functionality for other patch plotting options
% Figure out how to find outer approximation of 3D hybrid zonotope
% 
%	Trevor Bird - bird6@purdue.edu - Purdue 2021
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % 

function [ tCalcs ] = plot_MPT(obj,color,alpha,dims)
hold on;

% set defaults
if nargin == 1
	color = 'b';
	alpha = 1;
elseif nargin == 2
	alpha = 1;
end

Zh = obj;
Zh = getDimensions(Zh);
if Zh.n > 3 && nargin < 4
	error('Cannot plot in higher than 3 dimensions.. need to provide dims for 2D projection')
elseif Zh.n < 3 && nargin < 4
	dims = [ 1 2 ];
end

tic
% convert from HCG to H-rep
if Zh.ngb >= 1
	
	% if we have binaries we first need to find all feasible combinations
	if Zh.binTree
		% already did it online
		combos = Zh.binTree;
	else
		% solve MILP to find the binary tree
		combos = getLeaves(Zh);
	end
	tocingCombo = toc;
	
	% then use these combinations to convert to an equivalent complex of convex polyhedron
	% use MPT to convert to H-rep
	P = hyb2polyhedron_MPT(Zh,combos);
	tocingConvert = toc - tocingCombo;
	
else
	% no binaries, just find the single convex polyhedron
	tocingCombo = toc;

    P = hyb2polyhedron_MPT(Zh);

	tocingConvert = toc - tocingCombo;
	
end

% plot the complex of H-rep polyhedron

if P(1).Dim > 3
	P = P.projection([1,2]);
end

% find vertices and faces for using patch -> much faster than using MPT
for j = 1:numel(P)
    if size(P(j).V,1) ~= 1
        plot(P(j),'color',color,'alpha',alpha)
    end    
end

tocingPlot = toc - (tocingCombo + tocingConvert);

p = 1;

tocAll = toc;

tCalcs.combo = tocingCombo;
tCalcs.convert = tocingConvert;
tCalcs.plot = tocingPlot;
tCalcs.all = tocAll;

fprintf(['Plotted ',num2str(length(P),'%.0f'),' regions in ',num2str(tocAll,'%.2f'), ...
' sec\n  Finding feasible binaries: ',num2str((tocingCombo/tocAll)*100,'%.2f'),...
' percent \n  Converting to polyhedron: ',num2str((tocingConvert/tocAll)*100,'%.2f'),...
' percent \n  Plotting poly complex: ',num2str((tocingPlot/tocAll)*100,'%.2f'),' percent \n\n'])
	
end % end plot