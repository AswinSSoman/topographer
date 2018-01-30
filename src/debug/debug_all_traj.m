function draw_wgraph(graph, x2d, vertex_size, threshold, vertex_marker)
%DRAWTRANSITIONGRAPH Summary of this function goes here
%   Detailed explanation goes here

if ~exist('threshold', 'var'), threshold = 0; end
if ~exist('vertex_marker', 'var'), vertex_marker = {}; end

assert(iscell(vertex_marker));

% 1. draw the states (circle); do this before get axes informations
default_color = 'none';
for i = 1:numel(vertex_size);    
    s = vertex_size(i);
    rectangle('Position', [x2d(i, :) - s, [1 1] * s * 2], ...
        'Curvature', [1 1], ...
        'LineWidth', 2, ...
        'EdgeColor', [.7 .7 .7], ...
        'FaceColor', default_color);
end

hold on

% get coordinate width of the current axes
xlim = range(get(gca, 'XLim'));
ylim = range(get(gca, 'YLim'));
axes_size_x = norm([xlim, ylim]);

old_unit = get(gca, 'Units');

% get pixiel width of the current axes
set(gca, 'Units', 'pixels');
tmp = get(gca, 'Position');
axes_size_px = norm(tmp(3:4));

% get points width of the current axes
set(gca, 'Units', 'points');
tmp = get(gca, 'Position');
axes_size_pt = norm(tmp(3:4));

set(gca, 'Units', old_unit);  % restore the origin unit


% 2. draw the transitions (arrow)

% some drawing paramters
arrow_length = axes_size_px * 0.015;
arrow_width  = 1 / axes_size_pt * axes_size_px;  % about 2 points
bias_amount  = axes_size_x  * 0.006;
font_height = 10;  % dont change 
font_height_x = font_height / axes_size_pt * axes_size_x;
font_width_x = font_height_x * 0.5;
word_radius = font_width_x * 3.5 / 2; % 3 digits and 1 dot 

[from, to, weight] = find(graph);
for i = 1:numel(weight)
    f = from(i);
    t = to(i);
        
    % fot self-loop, put the weight inside the circle 
    if f == t
        
        marker = num2str(weight(i), '%.2f');
        if ~isempty(vertex_marker), 
            marker = [ '(' vertex_marker{f} ') ' marker ];
        else
        end
        
        text_x = x2d(f, :) - font_height_x * 0.1 - [word_radius 0];
        text(text_x(1), text_x(2), marker, ...
             'FontSize', font_height, 'FontUnit', 'points', ...
             'BackgroundColor', 'w', 'Margin', 2);
        continue;  % not doing any of the following 
    end
    
    % skip low weight links
    if weight(i) < threshold, continue; end
    
    f_x = x2d(f, :);
    t_x = x2d(t, :);
    v = x2d(t, :) - x2d(f, :);
    
    % trim the vector to avoid overlaping with the circle
    dirction = v / norm(v);
    f_x = f_x + dirction * vertex_size(f);
    t_x = t_x - dirction * vertex_size(t);

    
    % add a drift to the position, to avoid overlaping
    clockwise_rotate = [0 -1; 1, 0];
    bias_dir = dirction * clockwise_rotate';
    bias = bias_dir * bias_amount;
    
    % draw the arrow
    arrow(f_x + bias, t_x + bias, ...
        'Length', arrow_length, ...
        'BaseAngle', 45, ...
        'TipAngle', 22.5, ...
        'Width', arrow_width, ...
        'Color', [.7 .7 .7]);
    
    % annotate the arrow with edge weight
    if dirction(1) > 0, text_dir = - dirction;
    else text_dir = dirction; end
    text_x = (f_x + t_x)/2 - font_height_x * 0.1 + text_dir * word_radius;
    text_x = text_x + bias_dir * font_height_x / 2 + bias * 2;
    text(text_x(1), text_x(2), num2str(weight(i), '%.2f'), ...
         'FontSize', font_height, ...
         'FontUnit', 'points', ...
         'Rotation', atan(v(2)/v(1)) / pi * 180, ...
         'BackgroundColor', 'w', ...
         'Margin', 0.1);
end

end

