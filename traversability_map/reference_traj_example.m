clc; clear; close all;

% Define system parameters
dt = 0.1; % Time step
N = 10; % horizon length
T = N*dt; % Total time duration within horiozn

% Define initial and target states
x0 = [0; 0]; % current MPC state
x_target = [10; 10]; % from a straight line reference

% Generate reference trajectory using LQR
x_ref = generate_reference(x0, x_target, N, dt);
  
% Plot the reference trajectory
time = dt:dt:T;
figure;
subplot(2, 1, 1);
plot(time, x_ref(1, :), 'b', 'LineWidth', 2); hold on;
plot(time, x_ref(2, :), 'r', 'LineWidth', 2);
legend('x', 'y');
xlabel('Time [s]');
ylabel('States');
title('Reference State Trajectory');
grid on;