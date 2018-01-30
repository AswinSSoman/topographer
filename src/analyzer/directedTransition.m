function plot2DLandscape( x2d, adj, energy, ngrid, width, plot_more )
%PLOT2DLANDSCAPE Summary of this function goes here
%   Detailed explanation goes here

if ~exist('plot_more', 'var'), plot_more = false; end

x1 = x2d(:, 1);
x2 = x2d(:, 2);
padFactor = 0;
pad1 = range(x1) * padFactor; pad2 = range(x2) * padFactor;
x1lin = linspace(min(x1) - pad1, max(x1) + pad1, ngrid)';
x2lin = linspace(min(x2) - pad2, max(x2) + pad2, ngrid)';

% fill the blank with zeros energy points
newx1 = [];
newx2 = [];
x1width = x1lin(width) - x1lin(1);
x2width = x2lin(width) - x2lin(1);
for ix1 = x1lin(:)'
    x1l = ix1 - x1width;
    x1r = ix1 + x1width;
    
    for ix2 = x2lin(:)'
        x2l = ix2 - x2width;
        x2r = ix2 + x2width;
        
        if any( x1l < x1 & x1 < x1r & ...
                x2l < x2 & x2 < x2r )
            continue
        else
            newx1 = [newx1 ix1];
            newx2 = [newx2 ix2];
        end
        
    end
    
end
x1_filled = [x1; newx1'];
x2_filled = [x2; newx2'];
energy_filled = [energy(:); zeros(numel(newx1), 1)];

[X, Y] = meshgrid(x1lin, x2lin);
f = scatteredInterpolant(x1_filled, x2_filled, energy_filled, 'nearest');
Z = f(X, Y);

% make smoother; using a gaussian kernal
[H, V] = meshgrid((-width):width, (-width):width);
sigma = width/3;
weight = exp(- (H.^2 + V.^2) / (2 * sigma^2));
weight = weight / sum(weight(:));

[J, I] = meshgrid(1:ngrid, 1:ngrid);
Zpad = zeros(ngrid + 2*width);
Zpad(width+(1:ngrid), width+(1:ngrid)) = Z;
conv = @(i, j) sum(sum(weight .* Zpad(i:i+2*width, j:j+2*width)));
Z_smooth = arrayfun(conv, I, J);

% simplify
Z_smooth(Z_smooth >= -0) = inf;

hold on

% plot the mesh
colormap(brewermap(8, '*RdYlBu'));
%h = meshc(X, Y, Z_smooth);
surfc(X, Y, Z_smooth, 'EdgeColor', 'none', 'LineWidth', 1);

% plot cells and trajectory
if plot_more,
    % cells
    scatter3(x1, x2, energy(:), 25, [1 1 1] * 0.77, 'o');
    
    % trajectory
    [from, to, ~] = find(adj);
    for i = 1:numel(from)
        idx = [from(i) to(i)];
        plot3(x1(idx), x2(idx), [1 1] * -0.1, 'k', 'LineWidt', 1.5);
    end
end


hold off

% add styles
a = gca;
grid on
zlim([-1, 0]);
title('Landscape');
xlabel('x_1'); a.XTickLabel = [];
ylabel('x_2'); a.YTickLabel = [];
zlabel('Energy'); a.ZTick = [-1, -0.5, 0];
pbaspect([1 1 0.3]);
view(-30, 40);

end


function x = XTreeLayout(adj, pos)
root = 1;
offset = 0;
move = 0;
[idx, x] = getXTreeLayout(adj, pos, root, offset, move);
assert(numel(x) == size(adj, 1));

[~, ord] = sort(idx);
x = x(ord);
end

function [idx, x] = getXTreeLayout(adj, pos, root, offset, move)
next = find(any(adj(root, :), 1));
idx = root;
x = offset + move;
if ~isempty(next)
    
    % compute the offsets for subtree
    m = numel(next);
    if m == 1
        x_move = move * 0.8; % shrinking the move
    else
        x_move = pdist2(pos(next, :), pos(next(1), :));
        x_move = x_move - mean(x_move) + move * 0.7;
    end
    
    % position the subtree
    for i = 1:m
        [idx_i, x_i] = getXTreeLayout(adj, pos, next(i), offset+move, x_move(i));
        x = [x; x_i];
        idx = [idx; idx_i];
    end
end

end

