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