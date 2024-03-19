function b = CBF_sphere(states,obs)
% Inputs 
% states    : [x; y; z; psi; xdot; ydot; zdot; psidot] -> casadi variables;
% obs       : [obs_x; obs_y; obs_z; obs_r; obs_s] -> casadi variables
% outputs
% rho       : density function for spherical obstalce

    x_obs = obs(1);
    y_obs = obs(2);
    z_obs = obs(3);
    r_obs = obs(4);

    n_states = length(states);

    b = (states(1)-x_obs)^2 + (states(2)-y_obs)^2 + (states(3)-z_obs)^2 - r_obs^2;

end



