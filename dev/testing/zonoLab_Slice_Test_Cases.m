clear; clc;

%% Test 1
test_1 = zono([1 0 0 0; 0 2 0 0; 0 0 1 0; 0 0 0 3], zeros(4,1));
test_1_slice_A = zonoSlice(test_1, [0 0 1 0], 0);
test_1_slice_A = projection(test_1_slice_A, [1 2 4]);
figure; plot(test_1_slice_A, 'r', 0.2); grid on;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$x_4$','interpreter','latex');
title('Slice A $(x_3=0)$ of Test 1 Zono', 'interpreter', 'latex');

test_1_slice_B = zonoSlice(test_1, [0 0 0 1], 3);
test_1_slice_B = projection(test_1_slice_B, [1 2 3]);
figure; plot(test_1_slice_B, 'r', 0.2); grid on;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$x_3$','interpreter','latex');
title('Slice B $(x_4=3)$ of Test 1 Zono', 'interpreter', 'latex');
%% Test 2
test_2 = zono([1 0 1 0; 0 1 0 1; 0 0 1 1], [1; 0; -1]);
figure; plot(test_2, 'r', 0.2); grid on;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$x_3$','interpreter','latex');
title('Initial Test 2 Zono');

% test_2_slice_A = zonoSlice(test_2, [0 0 1], 0);
% figure; plot(test_2_slice_A, 'r', 0.2); grid on;
% xlabel('$x_1$','interpreter','latex');
% ylabel('$x_2$','interpreter','latex');
% title('Slice A of Test 2 Zono @ $x_3=0$', 'interpreter', 'latex');

test_2_slice_B = zonoSlice(test_2, [0 1 0], 0);
figure; [v,f] = plot(test_2_slice_B, 'r', 0.2); grid on;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$x_3$','interpreter','latex');
title('Slice B of Test 2 Zono @ $x_2=0$', 'interpreter', 'latex');

figure; [v_proj, f_proj] = plot(projection(test_2_slice_B, [1 3]), 'r', 0.2); grid on
%% Test 3
test_3 = conZono([1 0 1 0; 0 2 0 0.5; 0 0.5 1 0], zeros(3,1), [1 1 1 1], 0);
figure; [v3, f3] = plot(test_3, 'r', 0.2); grid on;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$x_3$','interpreter','latex');
title('Initial Test 3 ConZono');

test_3_slice = zonoSlice(test_3, [1 1 0], 0);
figure; [v3_slice, f3_slice] = plot(test_3_slice, 'r', 0.2); grid on;
xlabel('$x_1$','interpreter','latex');
ylabel('$x_2$','interpreter','latex');
zlabel('$x_3$','interpreter','latex');
title('Slice of Test 3 ConZono @ $x_1=-1$', 'interpreter', 'latex');

figure; plot(test_3_slice.projection(2:3), 'r', 0.2); grid on;
xlabel('$x_2$','interpreter','latex');
ylabel('$x_3$','interpreter','latex');
title('Projection of Test 3 Slice');