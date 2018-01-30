function group_transition = plotTransition(P, groups, x2d, rho)
%PLOTTRANSITION Summary of this function goes here
%   Detailed explanation goes here

axis equal
axis tight
hold on

assert(iscell(groups));
n = numel(groups);


% make group transition matrix
group_transition = zeros(n, n);
for i = 1:n
    for j = 1:n
        group_transition(i, j) = mean( sum( P(groups{i}, groups{j}) , 2 ) );
    end
end


% determinate width, depth and position of each group
group_width = zeros(n, 1);
group_depth = zeros(n, 1);
group_x     = zeros(n, 2);
for iter = 1:n 
    this_group = groups{iter};
    
    % position for this group 
    g_center_x = mean(x2d(this_group, :), 1);
    group_x(iter, :) = g_center_x; 

    % width of each group is the average distance to the center 
    group_width(iter) = mean( pdist2( x2d(this_group, :) , g_center_x ) ) / 2;
    
    % depth if defined by average density 
    group_depth(iter) = mean(rho(groups{iter}));
end


% % re-layout the vertex if time is provided  
% if exist('time', 'var'),     
%     x2d_new = zeros(size(x2d));
%     
%     [~, t_order] = sort(time(peaks), 'ascend');
%     head_group = t_order(1);
%     
%     forward_group = find(group_transition(head_group, :));
%     forward_group = setdiff(forward_group, head_group);
%     
% end


% plot the cells as reference 
for iter = 1:n 
    this_group = groups{iter};
    scatter(x2d(this_group, 1), x2d(this_group, 2), 36, 'filled');
end


% plot transition links 
dont_show_link = 0.001;   % filter out low weight links
draw_wgraph(group_transition, group_x, group_width, dont_show_link);

hold off
drawnow

end

