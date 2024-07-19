clc, close all, clear all

rng(1)

for i = 1:50
    Z{i} = randomSet(randi(2^32), 'hybZono', randi(10), randi(10), randi(10), randi(10));
    Z_rr{i} = removeRedundancy(Z{i});
end