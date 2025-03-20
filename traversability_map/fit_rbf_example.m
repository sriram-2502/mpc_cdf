%% fit rbfs for height maps
% Fit Gaussian RBFs to the height map data and plot the results using a fine mesh grid.
%
% Parameters:
% height_map - The height map (2D array) to fit.
% grid_spacing - The spacing between the centers in the uniform grid.
% sigma - The spread parameter for the Gaussian RBFs.
clc; clear; close all;

%% Example usage
grid_map = readmatrix('Terrain Map.xlsx','Sheet','center_hill');
max_val = max(grid_map,[],"all");
height_map = grid_map./max_val;

grid_spacing = 1;  % Adjust the spacing between RBF centers
sigma = 1;  % Tune sigma as needed

% Get the size of the height map
[rows, cols] = size(height_map);

% Generate the uniform grid centers
[center_x, center_y] = meshgrid(1:grid_spacing:cols, 1:grid_spacing:rows);
centers = [center_x(:), center_y(:)];

% Fit RBFs and get the weights
weights = fit_rbf(height_map, centers, sigma);

%%
% Generate the grid of coordinates
[x, y] = meshgrid(1:cols, 1:rows);
coordinates = [x(:), y(:)];

% Generate a fine grid of coordinates
[x_fine, y_fine] = meshgrid(linspace(1, cols, cols*10), linspace(1, rows, rows*10));
coordinates_fine = [x_fine(:), y_fine(:)];

% Reconstruct the height map using the RBFs and weights on the fine grid
num_centers = size(centers, 1);
rbf_matrix_fine = zeros(size(coordinates_fine, 1), num_centers);
for i = 1:num_centers
    diff = coordinates_fine - centers(i, :);
    rbf_matrix_fine(:, i) = exp(-sum(diff.^2, 2) / (2 * sigma^2));
end
fitted_height_map_fine = reshape(rbf_matrix_fine * weights, size(x_fine));

% Query a height value at point (x, y)
x_query = 11;  % Replace with your x-coordinate of the query point
y_query = 11;  % Replace with your y-coordinate of the query point
height_value = query_rbf_height([x_query, y_query], weights, centers, sigma);
% disp(['Height at (', num2str(x_query), ', ', num2str(y_query), ') is ', num2str(height_value)]);

%% Plot the fitted height map on a fine mesh grid
figure;
% Plot the original height map
subplot(1, 2, 1);
surf(x, y, height_map);
title('Original Height Map');
xlabel('X');
ylabel('Y');
zlabel('Height');
view(2)

subplot(1, 2, 2);
surf(x_fine, y_fine, fitted_height_map_fine, 'Edgecolor','none');
title('Fitted Height Map on Fine Mesh Grid');
xlabel('X');
ylabel('Y');
zlabel('Height');
% view(2)

% Make the figure look nicer
colormap jet;
colorbar;