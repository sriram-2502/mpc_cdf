function b = CBF_cylinder(states, obs)
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

    b = (sqrt((x(1)-cx)^2 + (x(2)-cz)^2))^2 - r^2;

end
