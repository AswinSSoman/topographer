function draw_density(X2d, rho)
%DRAWDENSITY Draw each cell, colored against its local density. Use this
%function as a building block for higher-level ploting. 
%
% @param X2d is the 2D coordinate for cells.
% @param rho is the local density for cells.

% set colors according to rho 
colormap = jet;
colornum = size(colormap, 1);
densityColor = colormap(max(ceil(colornum * rho/max(rho)), 1), :);

% plot all cells
scatter(X2d(:, 1), X2d(:, 2), 10, densityColor, 'filled');

end
