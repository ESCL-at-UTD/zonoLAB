% First, list all of the Clues that are given on the Sodoku board. If
% square (i,j) takes the value v, then append the triplet [i j v] to the
% matrix below. 

CLUES = [1 1 6
         1 2 5
         1 5 9
         2 1 4
         2 4 7
         2 7 3
         3 1 3
         3 6 1
         4 3 3
         4 6 6
         4 7 5
         4 8 4
         5 2 4
         5 5 5
         5 8 3
         6 2 9
         6 3 7
         6 7 6
         7 4 4
         7 9 2
         8 3 2
         8 6 7
         8 9 3
         9 5 1
         9 8 7
         9 9 5];

drawSudoku(CLUES)
title('Question')

%%

% For the most straight-forward first attempt, let's try the method of
% having 9 binary variables per square (this results in 729 binary 
% variables). This is inefficient, and could be reduced to as little as 4
% binary variables per square, but then the constraints become much more
% difficult to write. And I have found that Gurobi has no troubles handling
% the HZs that result from this code. 

Ac = [];
numCons = 4*81;
numHints = size(CLUES, 1);
Ab = zeros(numCons, 9^3);
thisCon = 0;

% the first numCons constraints come from the rules of solving a Sudoku
% puzzle, and the next numHints constraints are the hints given to us to
% begin the puzzle. 
b = [-7*ones(numCons, 1); ones(numHints,1)];

% each square gets a "sum-to-one" constraint (along the "value"-axis)
for i = 1:9
    for j = 1:9
        thisCon = thisCon+1;
        thisDims = [];
        for k = 1:9
            thisDims = [thisDims map(i,j,k)];
        end
        Ab(thisCon, thisDims) = 1;
    end
end

% each row gets a sum-to-one
for i = 1:9
    for k = 1:9
        thisCon = thisCon+1;
        thisDims = [];
        for j = 1:9
            thisDims = [thisDims map(i,j,k)];
        end
        Ab(thisCon, thisDims) = 1;
    end
end

% each column gets a sum-to-one
for j = 1:9
    for k = 1:9
        thisCon = thisCon+1;
        thisDims = [];
        for i = 1:9
            thisDims = [thisDims map(i,j,k)];
        end
        Ab(thisCon, thisDims) = 1;
    end
end

% each of the 3x3 mini squares gets a sum-to-one
for I = 1:3
    for J = 1:3
        for k = 1:9
            thisCon = thisCon+1;
            thisDims = [];
            for i = 1:3
                for j = 1:3
                    thisDims = [thisDims map(i+3*(I-1), j+3*(J-1), k)];
                end
            end
            Ab(thisCon, thisDims) = 1;
        end
    end
end

% and now here we constrain our solution to match the hints we were given
for m = 1:size(CLUES, 1)
    thisCon = thisCon+1;
    Ab(thisCon, map(CLUES(m,1), CLUES(m,2), CLUES(m,3))) = 1;
end

% time to make the hybzono; we don't actually care about the generators at
% all (we really just want to solve a feasibility problem), so here I just
% set them all zero. 
Gc = [];
Gb = zeros(1,9^3);
c = zeros(1,1);
Ac = zeros(numCons+numHints, 0);
Phi = hybZono(Gc, Gb, c, Ac, Ab, b);

%%

numLeaves = size(getLeaves(Phi,[]),2);
disp(['There is/are ' num2str(numLeaves) ' possible solution/s.'])

% if multiple leaves are found, then either: i) there is a bug in the code,
% or ii) the hints to the Sudoku board do not properly constrain you to a
% single solution. 
%%
SOL = [];
if numLeaves == 1
    XX = nan(9,9);
    xi = getLeaves(Phi);
    for i = 1:9
        for j = 1:9
            dims = [map(i,j,1) map(i,j,9)];
            val = find(xi(dims(1):dims(2))==1);
            SOL = [SOL; i j val];
        end
    end
    drawSudoku(SOL)
end
% you could edit this section to have it solve the sudoku board for all of
% the leaves (when there are multiple)
%%

% A triplet (i,j,k) corresponds with one possible value of one possible
% number on the Sudoku board. For example, (1,2,3) indicates the binary
% variable which determines whether or not square (1,2) takes the value 3.
% This function determines which of the 9^3 binary variables corresponds to
% the triplet (i,j,k)
function dim = map(i,j,k)
    dim = 81*(i-1) + 9*(j-1)+k;
end

% I found this function online
function drawSudoku(CLUES)
% Function for drawing the Sudoku board
%   Copyright 2014 The MathWorks, Inc. 

figure;hold on;axis off;axis equal % prepare to draw
rectangle('Position',[0 0 9 9],'LineWidth',3,'Clipping','off') % outside border
rectangle('Position',[3,0,3,9],'LineWidth',2) % heavy vertical lines
rectangle('Position',[0,3,9,3],'LineWidth',2) % heavy horizontal lines
rectangle('Position',[0,1,9,1],'LineWidth',1) % minor horizontal lines
rectangle('Position',[0,4,9,1],'LineWidth',1)
rectangle('Position',[0,7,9,1],'LineWidth',1)
rectangle('Position',[1,0,1,9],'LineWidth',1) % minor vertical lines
rectangle('Position',[4,0,1,9],'LineWidth',1)
rectangle('Position',[7,0,1,9],'LineWidth',1)

% Fill in the clues
%
% The rows of B are of the form (i,j,k) where i is the row counting from
% the top, j is the column, and k is the clue. To place the entries in the
% boxes, j is the horizontal distance, 10-i is the vertical distance, and
% we subtract 0.5 to center the clue in the box.
%
% If B is a 9-by-9 matrix, convert it to 3 columns first

if size(CLUES,2) == 9 % 9 columns
    [SM,SN] = meshgrid(1:9); % make i,j entries
    CLUES = [SN(:),SM(:),CLUES(:)]; % i,j,k rows
end

for ii = 1:size(CLUES,1)
    text(CLUES(ii,2)-0.5,9.5-CLUES(ii,1),num2str(CLUES(ii,3)))
end

hold off

end