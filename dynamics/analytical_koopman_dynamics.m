function [dx_dt] = analytical_koopman_dynamics(states, controls)
% Inputs
% states    : [1;x1; x2; x1^2] -> casadi variables;
% controls  : [u] -> casadi variables;
% Outputs
% dx_dt     : f(x)+g(x)u 
% F_sys     : f(x)
% G_sys     : g(x)

mu = 0.5;
lambda = 1.5;

x1 = states(1);
x2 = states(2);
u1 = controls(1);
u2 = controls(2);

% A_k = [1 0 0 0;0 mu 0 0;0 0 lambda -lambda;0 mu^2 0 0]; 
% B_k = [0;0;1;0];
% A_k_reduced = [mu 0 0;0 lambda -lambda;mu^2 0 0]; 
% B_k_reduced = [0;1;0];
% dx_dt = A_k_reduced*[z2;z3;z4] + B_k_reduced*u;

x1_dot = mu*x1 + u1;
x2_dot = lambda*(x2*x1^2) + u2;

dx_dt = [x1_dot;x2_dot];