function b = CBF_circle(states,obs)
% Inputs 
% states    : [x; y; z; psi; xdot; ydot; zdot; psidot] -> casadi variables;
% obs       : [obs_x; obs_y; obs_z; obs_r] -> casadi variables
% outputs
% b(x)       : cbf function for circular obstalce

    x_obs = obs(1);
    y_obs = obs(2);
    r_obs = obs(3);
    b = (states(1)-x_obs)^2 + (states(2)-y_obs)^2 - r_obs^2;

end



