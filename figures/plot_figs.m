% List figures (including 3D and 2D plots)
figureList3D = {'example1_cdf.fig', 'example2_cbf_vs_cdf.fig'};
figureList2D = {'example2_cbf_vs_cdf_logs.fig'};

% Total number of subplots
totalSubplots = numel(figureList3D) + numel(figureList2D);

% List new line colors for each figure
lineCol3D = jet(numel(figureList3D));
lineCol2D = jet(numel(figureList2D));

% Create a new figure
fh = figure();

% Plot 3D figures
for i = 1:numel(figureList3D)
    % Open figure-i
    fighand = openfig(figureList3D{i}, 'invisible'); % Open invisibly
    % Get axis handle (assumption: only 1 axes in figure)
    axHand = findobj(fighand, 'Type', 'Axes');
    
    % Create a subplot for this figure
    newAxes = subplot(2, 4, i, 'Parent', fh); % Adjust based on totalSubplots
    hold(newAxes, 'on');
    grid(newAxes, 'on');
    
    % Copy content from axHand to newAxes
    for obj = axHand.Children'
        h = copyobj(obj, newAxes);
        % Set line color
        % if isprop(h, 'Color')
        %     h.Color = lineCol3D(i, :); % Set the line color
        % end
    end
    
    % Set axis properties for 3D plots
    colormap gray
    view(newAxes, [-23, 16]); % Set the 3D view
    box(newAxes, 'on');       % Add a box around the plot
    axis(newAxes, 'square');  % Square axis
    set(newAxes, 'FontSize', 15, 'LineWidth', 1.5); % Font size and line width
    
    % Set axis labels with LaTeX formatting
    xlabel(newAxes, '$x_1$', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(newAxes, '$x_2$', 'Interpreter', 'latex', 'FontSize', 20);
    zlabel(newAxes, '$x_3$', 'Interpreter', 'latex', 'FontSize', 20);
    
    % Set the title of the subplot
    title(newAxes, axHand.Title.String, 'Interpreter', 'None');
    
    % Close original figure
    close(fighand);
end

% Plot 2D figures
for i = 1:numel(figureList2D)
    % Open figure-i
    fighand = openfig(figureList2D{i}, 'invisible'); % Open invisibly
    % Get axis handle (assumption: only 1 axes in figure)
    axHand = findobj(fighand, 'Type', 'Axes');
    
    % Create a subplot for this figure
    newAxes = subplot(2, 4, [5, 6], 'Parent', fh); % Adjust based on totalSubplots
    hold(newAxes, 'on');
    grid(newAxes, 'on');
    
    % Copy content from axHand to newAxes
    for obj = axHand.Children'
        h = copyobj(obj, newAxes);
        % % Set line color
        % if isprop(h, 'Color')
        %     h.Color = lineCol2D(i, :); % Set the line color
        % end
    end
    
    % Set axis properties for 2D plots
    box(newAxes, 'on');       % Add a box around the plot
    set(newAxes, 'FontSize', 15, 'LineWidth', 1.5); % Font size and line width
    
    % Set axis labels with LaTeX formatting
    xlabel(newAxes, 'time, (s)', 'Interpreter', 'latex', 'FontSize', 20);
    ylabel(newAxes, 'distance', 'Interpreter', 'latex', 'FontSize', 20);
    

    % Set the title of the subplot
    title(newAxes, axHand.Title.String, 'Interpreter', 'None');
    
    % Close original figure
    close(fighand);
end
