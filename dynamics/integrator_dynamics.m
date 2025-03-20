function [dx_dt,F_sys,G_sys] = integrator_dynamics(states, controls)
% Inputs
% states    : [x; y; theta; v] -> casadi variables;
% controls  : [v;a] -> casadi variables;
% Outputs
% dx_dt     : f(x)+g(x)u 
% F_sys     : f(x)
% G_sys     : g(x)

x = states(1);
y = states(2);

u1 = controls(1);
u2 = controls(2);

F_sys = [0;
         0];

G_sys = [1 0;
         0 1];

dx_dt = [u1;
         u2];