function rho = density_cylinder(states, obs)
% Inputs 
% states    : [x; y; z; psi; xdot; ydot; zdot; psidot] -> casadi variables;
% obs       : [obs_x; obs_y; obs_z; obs_r; obs_s] -> casadi variables
% outputs
% rho       : density function for spherical obstalce

    cx = obs(1); cy = obs(2); cz = obs(3); % center
    r = obs(4); % radius of cylinder
    r1 = r; r2 = obs(5); % for sensing

    n_states = length(states);
    P_lyap = eye(n_states);
    alpha = 0.1;

    shape = (sqrt((x(1)-cx)^2 + (x(2)-cz)^2))^2 - r^2;
    shift = 1; % to get rid of random zeros
    h =  shape + shift;
    temp1 = h/(r2^2-r1^2);

    f_bar1 = if_else(temp1 > 0, exp(-1/temp1), 0)/(if_else(temp1 > 0, exp(-1/temp1), 0) + if_else(1-temp1 > 0, exp(-1/(1-temp1)), 0));
    Phi = f_bar1;
    V = [states(1) states(2) states(3) states(4) states(5) states(6) states(7) states(8)] *P_lyap*([states(1); states(2); states(3); states(4); states(5); states(6); states(7); states(8)]);
    rho = Phi/(V.^alpha);
end
