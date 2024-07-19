clc; clear; close all;
%% needs work. does not meet requirements!!


% Define the polynomial reference trajectory as an inline function
a = 0.5;
b = 0;
c = 0;
d = 0;
y_ref_func = @(x) a * x.^3 + b * x.^2 + c * x + d;

% Define the range of values for the current state
x_range = -3:0.01:3;
y_range = -3:0.01:3;

% Initialize a matrix to store the penalty values
penalty_values = zeros(length(x_range), length(y_range));

% Define epsilon and k
epsilon = 1;
k = 100;

% Define a point to query
x_query = 1.0;
y_query = y_ref_func(x_query) + 0.3;  % Within epsilon distance

% Calculate the y value of the polynomial at x_query using the inline function
y_ref_query = y_ref_func(x_query);

% Calculate the distance to the polynomial curve
distance_query = abs(y_query - y_ref_query);

% Calculate the penalty based on the distance
penalty_query = penaltyFunction(distance_query, epsilon, k);

% Display the results
fprintf('Query Point: (%.2f, %.2f)\n', x_query, y_query);
fprintf('Distance to Reference: %.2f\n', distance_query);
fprintf('Penalty: %.2f\n', penalty_query);

% Calculate the penalty for each value in the range
for i = 1:length(x_range)
    for j = 1:length(y_range)
        x = x_range(i);
        y = y_range(j);

        % Calculate the y value of the polynomial at x
        a=1; b=2;
        c=1; d=1;
        y_ref = y_ref_func(x);

        % Calculate the distance to the polynomial curve
        distance = abs(y - y_ref);

        % Calculate the penalty based on the distance
        penalty_values(j, i) = -penaltyFunction(distance, epsilon, k);
    end
end
penalty_values = penalty_values - min(penalty_values);
% Plot the penalty function
[X, Y] = meshgrid(x_range, y_range);
surf(X, Y, penalty_values, 'EdgeColor', 'None');
hold on;

% Plot the reference polynomial trajectory
x_ref_range = -2:0.1:2;
y_ref_range = y_ref_func(x_ref_range);
plot3(x_ref_range, y_ref_range, zeros(size(x_ref_range)), 'r-', 'LineWidth', 2);

xlabel('x');
ylabel('y');
zlabel('Penalty');
title('Penalty Function for Polynomial Reference Trajectory');
legend('Penalty Surface', 'Reference Trajectory');
hold off;

%% tracking function
function penalty = penaltyFunction(distance, epsilon, k)
    % Penalty function that takes zero value within a tube around the reference trajectory
    % and becomes large outside.
    %
    % Parameters:
    % distance (scalar): Distance from the current state to the reference trajectory
    % epsilon (scalar): Width of the tube around the reference trajectory
    % k (scalar): Scaling factor for the penalty outside the tube
    %
    % Returns:
    % penalty (scalar): Penalty value

    % Check if the distance is within the specified tube
    if distance >= epsilon
        penalty = 0;
    else
        penalty = k * (distance - epsilon)^2;
    end
end