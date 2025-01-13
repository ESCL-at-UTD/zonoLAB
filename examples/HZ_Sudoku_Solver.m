% for the first attempt, let's try the method of just 9 binary variables
% per square (this results in 729 binary variables, which is bad)

% each square have binary variables {xi_(i,j,1),...,xi_(i,j,9)}

% also, I'm going to write all the hybZonos in \xi_b \in {0, 1} form, and
% then convert it to the usual {-1, 1} later.

Ac = [];
numCons = 4*81;
Ab = zeros(numCons, 9^3);
thisCon = 0;

b = ones(numCons, 1);

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

Gc = [];
Gb = [];
c = zeros(1,1);

Phi = hybZono([], Gb, c, [], Ab, b);




%%

function dim = map(i,j,k)
    dim = 81*(i-1) + 9*(j-1)+k;
end