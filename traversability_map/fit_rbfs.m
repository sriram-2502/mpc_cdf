%% fit rbfs for height maps
clc; clear; close all;

%% Example usage
grid_map = readmatrix('Terrain Map.xlsx','Sheet','C_shaped');
max_val = max(grid_map,[],"all");
height_map = grid_map./max_val;

grid_spacing = 3;  % Adjust the spacing between RBF centers
sigma = 1;  % Replace with your desired sigma
rbf_fitting_and_plot_fine_mesh(height_map, grid_spacing, sigma);

%% plotting functions
function rbf_fitting_and_plot_fine_mesh(height_map, grid_spacing, sigma)
    % Fit Gaussian RBFs to the height map data and plot the results using a fine mesh grid.
    %
    % Parameters:
    % height_map - The height map (2D array) to fit.
    % grid_spacing - The spacing between the centers in the uniform grid.
    % sigma - The spread parameter for the Gaussian RBFs.

    % Get the size of the height map
    [rows, cols] = size(height_map);
    
    % Generate the uniform grid centers
    [center_x, center_y] = meshgrid(1:grid_spacing:cols, 1:grid_spacing:rows);
    centers = [center_x(:), center_y(:)];
    
    % Fit RBFs and get the weights
    weights = fit_rbf(height_map, centers, sigma);
    
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
    
    % Plot the fitted height map on a fine mesh grid
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
    surf(x_fine, y_fine, fitted_height_map_fine);
    title('Fitted Height Map on Fine Mesh Grid');
    xlabel('X');
    ylabel('Y');
    zlabel('Height');
    view(2)
    
    % Make the figure look nicer
    colormap jet;
    colorbar;
end

%% fit rbf function 
function weights = fit_rbf(height_map, centers, sigma)
    % Fit Gaussian RBFs to the height map data.
    % 
    % Parameters:
    % height_map - The height map (2D array) to fit.
    % centers - The coordinates of the centers of the RBFs (2D array of shape [num_centers, 2]).
    % sigma - The spread parameter for the Gaussian RBFs.
    %
    % Returns:
    % weights - The weights for the RBFs.

    % Get the size of the height map
    [rows, cols] = size(height_map);
    
    % Generate the grid of coordinates
    [x, y] = meshgrid(1:cols, 1:rows);
    coordinates = [x(:), y(:)];
    
    % Number of centers
    num_centers = size(centers, 1);
    
    % Compute the RBF matrix
    rbf_matrix = zeros(rows * cols, num_centers);
    for i = 1:num_centers
        diff = coordinates - centers(i, :);
        rbf_matrix(:, i) = exp(-sum(diff.^2, 2) / (2 * sigma^2));
    end
    
    % Flatten the height map
    heights = height_map(:);
    
    % Solve the linear system to find the weights
    weights = (rbf_matrix' * rbf_matrix) \ (rbf_matrix' * heights);
end