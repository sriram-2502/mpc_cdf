function b = CBF_torous_xz(states)
% Inputs 
% states    : [x; y; z; psi; xdot; ydot; zdot; psidot] -> casadi variables;
% outputs
% rho       : density function for spherical obstalce
    
    % torus 1 in XZ
    cx = obs(1); cy = obs(2); cz = obs(3); % center
    R = 7; % distance from center
    r = obs(4); % radius of cylinder
    r1 = r; r2 = obs(5); % for sensing

    n_states = length(states);
    P_lyap = eye(n_states);
    alpha = 0.1;

    b = (sqrt((states(1)-cx)^2 + (states(2)-cy)^2) -R)^2 + (states(3)-cz)^2 - r^2;
end
