function plotGeneProfiles( time, expr, gene_ids, ltree)
% plotGeneProfiles plot a smoothed expression profile for each row of 'expr'
% matrix. 'ltree' must be provided to determinate the branching structure of
% all cells; the branching is visualized as different style of lines. Also,
% each gene profile is drawn in different color. 
%
% INPUT:
%   time: the psedotime for each cell (coloum in 'expr' matrix).
%   expr: an N-by-M expression data matrix for N cells and M genes. 
%   gene_ids: labels for each of n genes. 
%   ltree: a tree-structure to represent the branching relationship of cells; 
%     this data is computed by linkedTree() function. 
%

% verify the inpyt format 
assert(numel(time) == size(expr, 1));
assert(size(expr, 2) == numel(gene_ids));

% get all branches by processing the LTree. 
subname_f = @(i) [ 'subtree' num2str(i) ];
all_branches = split_tree(ltree, subname_f, 'members');


%{
% % use 'plot_as_function' for simple Y-branching 
% if numel(all_branches) == 2
%     branch1 = all_branches{1};
%     branch2 = all_branches{2};
%     trunk = intersect(branch1, branch2);
%     
%     % marking in the fashion of Wanderlust 
%     branchY = zeros(1, numel(time));
%     branchY(branch1) = 1;
%     branchY(branch2) = -1;
%     branchY(trunk) = 0;
%     
%     plot_as_function(time, expr, 'branchY', branchY, 'labels', gene_ids, ...
%         'show_error', true, 'smooth', 0.1);
%     
%     return ;
% end
%}


%%%
% For more complex, multiple branching, we develop a similar method. 

num_cells = numel(time);
num_genes = numel(gene_ids);

% normalize the expr data 
for i = 1:num_genes
    expr(:, i) = rescale01(expr(:, i));
end

% normalize t, and define bin time
time = rescale01(time); % normalize to 0~1
bin_num = 100;  % the number of points in the final smoothed plot 
bin_time = linspace(0, 1, bin_num);


%%%
% global smoothing matrix; initially taking weights from all branches for each
% bin. But we would modify it for each branch, individually. 
smooth_mat = zeros(bin_num, num_cells);
for i = 1:bin_num  % each row is a kernel 
    smooth_mat(i, :) = smooth_weight(bin_time(i), time, 0.05);
end


% spercific the styling of plot lines 
gene_color = brewermap(num_genes, 'Set1'); % color from Brewermap package 
branch_style = {'-', ':', '--', '-.'};
if numel(all_branches) > numel(branch_style)
    disp('Too much branches; the plot might be confusing');
end

% plot all genes for each branch; also modified the smoot_mat. 
decay_factor = 0.0; % how long the influence from other branches persist
for b = 1:numel(all_branches)
    %%%
    % Eliminating weights from other branches; but we do it only partially, in
    % order to achieve a smooth transition from branches. 
    smooth_mat_modified = nan(size(smooth_mat));
    weight_limit = min( smooth_mat(:, all_branches{b}), [], 2);
    for k = 1:numel(all_branches)
        if k == b, 
            % just copy the weights; in-branch weights dont need modification 
            smooth_mat_modified(:, all_branches{b}) = smooth_mat(:, all_branches{b});
            continue; 
        end
        
        % decay the weights of cells from other brancehes 
        different_cells = setdiff(all_branches{k}, all_branches{b});
        decay = rescale01(time(different_cells)) .^ decay_factor;
        reweight = repmat(1 - decay, bin_num, 1);
        smooth_mat_modified(:, different_cells) = ...
            min( smooth_mat(:, different_cells) .* reweight, ...
                 repmat(weight_limit, 1, numel(different_cells)) );
    end
    
    % normalizing factor for the smooth_mat weights 
    norm_factor = repmat( sum(smooth_mat_modified, 2) , 1, numel(gene_ids));
    normed_expr = ( smooth_mat_modified * expr ) ./ norm_factor;
    
    % TODO: plot the only the new part
    new_cells = [];
    
    hold on; % must placed before set(, "ColorOrder") stuffs
    set(gca, 'ColorOrder', gene_color);
    set(gca, 'ColorOrderIndex', 1);
    plot(bin_time, normed_expr, ...
        'LineWidth', 3, ...
        'LineStyle', branch_style{b});    
end


% add legend and other styling-stuffs
if exist('gene_ids', 'var')
    legend(gene_ids, 'FontSize', 10);
end
xlabel('Pseudotime'); 
ylabel('Expression');
title('Gene Profiles');
xlim([0 1]);
set(gca, 'FontSize', 14);

hold off; 
drawnow;

end

function w = smooth_weight(x0, x, factor)
% use gaussian kernal for smoothing
%
% w is un-normalized ! 
%
% TODO: more smoothing options 

w =  exp(- ((x - x0)/factor) .^2  ); 

end

function x = rescale01(x)
if isempty(x), return; end
x = ( x - min(x) ) / range(x);
if range(x) == 0, error('rescale01: range 0.'); end
end
