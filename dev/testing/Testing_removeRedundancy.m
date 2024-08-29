clc, close all, clear all

flag_plots = false;
which_to_plot = [];
rng(2)

data = nan(0,3);

num_trials = 20;
attempts = nan(num_trials, 1);
Z = repmat({nan}, num_trials, 1);
Z_rr = Z; 
times = nan(num_trials, 1);
for i = 1:num_trials
    Z{i} = randomSet(randi(2^32), 'hybZono', 10, 50, 10, 30);
    attempts(i) = 1;
    while is_empty(Z{i})
        Z{i} = randomSet(randi(2^32), 'hybZono', 10, 50, 10, 30);
        attempts(i) = attempts(i)+1;
    end
%     Z{i}.b = [Z{i}.b; 2*Z{i}.b];
%     Z{i}.Ab = [Z{i}.Ab; 2*Z{i}.Ab];
%     Z{i}.Ac = [Z{i}.Ac; 2*Z{i}.Ac];
    tic
    Z_rr{i} = removeRedundancy(Z{i});
    times(i) = toc;
    M = eye(3, Z{i}.n);
    disp(i)
    
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

table(data(:,1), data(:,2), data(:,3), times, 'VariableNames', {'i', 'nC (orig)', 'nC (red)', 'seconds'})

return
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
hold on, grid on, grid minor, axis equal
plot(Z, 'r', 0.6)

%%
clc, close all, clear all

Z = randomSet(1, 'hybZono', 10, 20, 5, 4);
Z.Ac = [Z.Ac; zeros(1, Z.nGc)];
Z.Ab = [Z.Ab; zeros(1, Z.nGb)];
Z.b = [Z.b; 0];

Ac_new = eye(1, Z.nGc);
Ab_new = zeros(1, Z.nGb);
b_new = 0.5;

Z.Ac = [Z.Ac(1:2,:); Ac_new; Z.Ac(3:5,:)];
Z.Ab = [Z.Ab(1:2,:); Ab_new; Z.Ab(3:5,:)];
Z.b = [Z.b(1:2,:); b_new; Z.b(3:5,:)];

Z_simp = removeRedundancy(Z);