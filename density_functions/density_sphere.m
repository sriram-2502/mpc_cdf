function rho = density_sphere(states,obs)
% Inputs 
% states    : [x; y; z; psi; xdot; ydot; zdot; psidot] -> casadi variables;
% obs       : [obs_x; obs_y; obs_z; obs_r; obs_s] -> casadi variables
% outputs
% rho       : density function for spherical obstalce

    x_obs = obs(1);
    y_obs = obs(2);
    z_obs = obs(3);
    r_obs = obs(4);
    s_obs = obs(5);

    n_states = length(states);
    P_lyap = eye(n_states);
    alpha = 0.1;
    
    shape = (states(1)-x_obs)^2 + (states(2)-y_obs)^2 + (states(3)-z_obs)^2 - r_obs^2;
    temp1 = shape/(s_obs^2-r_obs^2);

    f_bar1 = if_else(temp1 > 0, exp(-1/temp1), 0)/(if_else(temp1 > 0, exp(-1/temp1), 0) + if_else(1-temp1 > 0, exp(-1/(1-temp1)), 0));
    Phi = f_bar1;
    V = [states(1) states(2) states(3) states(4) states(5) states(6) states(7) states(8)] *P_lyap*([states(1); states(2); states(3); states(4); states(5); states(6); states(7); states(8)]);
    rho = Phi/(V.^alpha);
end



