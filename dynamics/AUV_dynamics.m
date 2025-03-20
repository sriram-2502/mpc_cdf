function [dx_dt,F_sys,G_sys] = AUV_dynamics(states, controls, dt)
% Inputs
% states    : [x; y; z; psi; xdot; ydot; zdot; psidot] -> casadi variables;
% controls  : [u;v;w;r] -> casadi variables;
% Outputs
% dx_dt     : f(x)+g(x)u 
% F_sys     : f(x)
% G_sys     : g(x)

% ------ unpack states and control -------------
x = states(1); xdot = states(5);   
y = states(2); ydot = states(6);
z = states(3); zdot = states(7);
psi = states(4); psidot = states(8); 

u = controls(1); v = controls(2);
w = controls(3); r = controls(4);

 %---- Parameters setup -----------------
        m = 54.54; %Kg
        G = 535; %N
        B = 53.4;
        Iz = 13.587;
        Xuu = 2.3e-2;
        Yvv = 5.3e-2;
        Zww = 1.7e-1;
        Nrr = 2.9e-3;
        Xudot = -7.6e-3;
        Yvdot = -5.5e-2;
        Zwdot = -2.4e-1;
        Nrdot = -3.4e-3;
        Xu = 2e-3;
        Yv = -1e-1;
        Zw = -3e-1;
        Nr = 3e-2;

       J = [cos(psi) -sin(psi) 0 0;
            sin(psi)   cos(psi) 0 0;
            0 0 1 0;
            0 0 0 1];
       X1_dot = [xdot; ydot; zdot; psidot];

       g = [0; 0; -(G-B); 0];
       Jdot = [-psidot*sin(psi) -psidot*cos(psi) 0 0;
                psidot*cos(psi) -psidot*sin(psi) 0 0;
                0 0 0 0;
                0 0 0 0];
       invJ = [cos(psi) sin(psi) 0 0;
               -sin(psi) cos(psi) 0 0;
               0 0 1 0;
               0 0 0 1];
       m11 = m - Xudot;
       m22 = m - Yvdot;
       m33 = m - Zwdot;
       m44 = Iz - Nrdot;
       Mt = diag([m11 m22 m33 m44]);
       M = (invJ)' * Mt * invJ;
       invM = [(7676892219811089*cos(psi)*cos(psi))/140737488355328+(10919*sin(psi)*sin(psi))/200,  (7676892219811089*cos(psi)*sin(psi))/140737488355328-(10919*sin(psi)*cos(psi))/200,       0,        0;
                    (7676892219811089*sin(psi)*cos(psi))/140737488355328-(10919*cos(psi)*sin(psi))/200,  (10919*cos(psi)*cos(psi))/200+(7676892219811089*sin(psi)*sin(psi))/140737488355328,       0,        0;
                                                                                                           0,                                                                                        0, 2739/50,        0;
                                             0,                                                                                        0,       0, 8494/625];                                                                                           
        
       V = invJ * [xdot; ydot; zdot; psidot];
        d11 = -Xu - Xuu * abs(V(1));
        d22 = -Yv - Yvv * abs(V(2));
        d33 = -Zw - Zww * abs(V(3));
        d44 = -Nr - Nrr * abs(V(4));                                                                                
        Dv = diag([d11, d22, d33, d44]);
        D = (invJ)' * Dv * invJ;
        
        
        Cv = [0 0 0 -(m-Yvdot)*V(2);
                  0 0 0 (m-Xudot)*V(1);
                  0 0 0 0;
                  (m-Yvdot)*V(2) -(m-Xudot)*V(1) 0 0];
        C = (invJ)' * (Cv - M * invJ * Jdot) * invJ;
        fx = invM * (-C * [xdot; ydot; zdot; psidot] - D * [xdot; ydot; zdot; psidot]);
        F_sys = [xdot; ydot; zdot; psidot; fx];   
        G_sys = dt*[zeros(4);invM * (invJ)'];
        
        u_sys = invM * (invJ)' * controls;
        X2_dot = fx + u_sys;
        
    dx_dt = [X1_dot; X2_dot];
end