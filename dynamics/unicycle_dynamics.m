function [dx_dt,F_sys,G_sys] = unicycle_dynamics(states, controls, mismatch)
% Inputs
% states    : [x; y; theta; v] -> casadi variables;
% controls  : [v;a] -> casadi variables;
% Outputs
% dx_dt     : f(x)+g(x)u 
% F_sys     : f(x)
% G_sys     : g(x)

x = states(1);
y = states(2);
theta=states(3);
v = states(4);

w = controls(1);
a = controls(2);

if(mismatch)
    m = 0.1;
    a = controls(2)/m;
end

F_sys = [v.* cos(theta);
         v.*sin(theta);
         0;
         0];

G_sys = [0 0;
         0 0;
         1 0;
         0 1];

dx_dt = [v.* cos(theta);
         v.*sin(theta);
         w;
         a];

    