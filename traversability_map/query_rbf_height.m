function height = query_rbf_height(point, weights, centers, sigma)
    % Query the height value at a given point (x, y) using the fitted RBFs.
    %
    % Parameters:
    % point - The [x,y]-coordinate of the query point (1,2).
    %
    % Returns:
    % height - The height value at the query point.
    import casadi.*
    num_centers = size(centers, 1);
    rbf_values = casadi.SX.zeros(1, num_centers);
    for i = 1:num_centers
        diff = point - centers(i, :);
        rbf_values(i) = exp(-sum(diff.^2, 2) / (2 * sigma^2));
    end
    
    height = rbf_values * weights;
end