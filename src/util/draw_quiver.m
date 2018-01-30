function draw_quiver(from_x, from_y, to_x, to_y, size, draw_head, color)
% a little function for drawing an arrow.

if ~exist('size', 'var'), size = 10; end
if ~exist('draw_head', 'var'), draw_head = true; end
if ~exist('color', 'var'), color = 'k'; end

if draw_head
    head_option = 'on';
else
    head_option = 'off';
end

quiver(from_x, from_y, to_x - from_x, to_y - from_y, 0,...
    'LineWidth', size, ...
    'MaxHeadSize', size * 10, ...
    'ShowArrowHead', head_option, ...
    'Color', color);

end