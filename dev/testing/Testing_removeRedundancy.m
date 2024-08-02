clc, close all, clear all

flag_plots = false;
which_to_plot = [8];
rng(2)

data = nan(0,3);

num_trials = 10;
attempts = nan(num_trials, 1);
Z = repmat({nan}, num_trials, 1);
Z_rr = Z;
for i = 1:num_trials
    Z{i} = randomSet(randi(2^32), 'hybZono', 5, randi(10), randi(10), randi(10));
    attempts(i) = 1;
    while is_empty(Z{i})
        Z{i} = randomSet(randi(2^32), 'hybZono', 5, randi(10), randi(10), randi(10));
        attempts(i) = attempts(i)+1;
    end
    Z_rr{i} = removeRedundancy(Z{i});
    M = eye(3, Z{i}.n);
    
    if flag_plots & any(i==which_to_plot)
        figure
        subplot 121
            hold on, grid on, grid minor
            title(['Random HZ (' num2str(i) ')'])
            subtitle(['$n_C = ' num2str(Z{i}.nC) '$'], 'interpreter', 'latex')
            plot(M*Z{i}, 'b', 0.6)

        subplot 122
            hold on, grid on, grid minor
            title 'Reduced HZ'
            subtitle(['$n_C = ' num2str(Z_rr{i}.nC) '$'], 'interpreter', 'latex')
            plot(M*Z_rr{i}, 'g', 0.6)
    end
    
    data(i,:) = [i Z{i}.nC Z_rr{i}.nC];
end

table(data(:,1), data(:,2), data(:,3), 'VariableNames', {'i', 'nC (orig)', 'nC (red)'})
%%

clc, close all, clear all

Gc = [1 2 0;
      0 0 0;
      0 0 1];
Gb = [0;0;0];
c = zeros(3,1);
Ac = [0 1 0];
Ab = 0;
b = [0];
Z = hybZono(Gc, [], c, Ac, [], b);

% Z = hybZono(Gc, [], c, [], [], []);
figure
plot(Z, 'r', 0.6)
hold on, grid on, grid minor, axis equal