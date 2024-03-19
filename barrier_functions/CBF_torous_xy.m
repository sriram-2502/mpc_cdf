function b = CBF_torous_xy(states, obs)
% Inputs 
% states    : [x; y; z; psi; xdot; ydot; zdot; psidot] -> casadi variables;
% obs       : [obs_x; obs_y; obs_z; obs_r] -> casadi variables
% outputs
% b(x)       : cbf function for toroz in xy plane
    
    % torus 1 in XY
    cx = obs(1); cy = obs(2); cz = obs(3); % center
    R = 10; % distance from center
    r = obs(4); % radius of cylinder
    b = (sqrt((states(1)-cx)^2 + (states(2)-cy)^2) -R)^2 + (states(3)-cz)^2 - r^2;

end

