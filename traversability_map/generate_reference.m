function x_ref = generate_reference(x0, x_target, N, dt)
    % Generates a reference trajectory using LQR control for a 2D single integrator
    %
    % Inputs:
    % x0 - Initial state [x; y]
    % x_target - Target state [x; y]
    % T - Total time duration
    % dt - Time step
    %
    % Outputs:
    % time - Time vector
    % x_ref - Reference state trajectory
    % u_ref - Reference control input trajectory

    % Define the system matrices for the 2D single integrator
    A = eye(2);
    B = dt * eye(2);
    
    % Define LQR weight matrices
    Q = diag([100, 100]); % State weighting matrix
    R = diag([1, 1]); % Input weighting matrix

    % LQR design
    [K, ~, ~] = dlqr(A, B, Q, R);

    % Initialize state and input trajectories
    num_steps = N;
    nx = size(A, 1);
    nu = size(B, 2);
    x_ref = zeros(nx, num_steps);
    u_ref = zeros(nu, num_steps - 1);

    % Set initial state
    x_ref(:, 1) = x0;

    % Simulate the system using the LQR controller
    for k = 1:num_steps - 1
        % Calculate the control input
        u_ref(:, k) = -K * (x_ref(:, k) - x_target);

        % Update the state
        x_ref(:, k + 1) = A * x_ref(:, k) + B * u_ref(:, k);
    end
end