function draw_branching(X2d, progen, trajectory_mat)
%DRAWBRANCHING Draw a colored branching tree on the 2d cell space. Use this
%function as a building block for higher-level ploting


% apply gray-diant on the cell color, for better resolution
%cellColor = cellColorMixed .* repmat(0.2 + pdt(:) * 0.8, 1, 3);

% draw the trajectory, as a sequence of arrows
hold on
markersize = 4;
draw_branch_recursive(X2d, progen, trajectory_mat, markersize);

% mark each branch end
% branchArray = cell(1, nBranches);
% for i = 1:nBranches
%     branchArray{i} = S.(['branch' num2str(i)]);
% end
% branchnames = cell(nBranches, 1);
% for i = 1:nBranches
%     branchnames{i} = ['branch' num2str(i)];
%     branchEnd = branchArray{i}.cellIndices(end);
%     scatter(X2d(branchEnd, 1), X2d(branchEnd, 2), 25, ...
%         branchColor(i, :), 'filled');
% end
% legend(branchnames, 'FontSize', 16, 'FontWeight', 'bold');


end % end function

function branch_num = draw_branch_recursive(x, root, adj_mat, msize)

leaves = find(adj_mat(root, :));
branch_num = numel(leaves);

if isempty(leaves), return; end

for leaf = leaves(:)'
    
    % draw for leaf nodes first
    sub_branch_num = draw_branch_recursive(x, leaf, adj_mat, msize);
    
    if sub_branch_num == 1 % intermediate segments, no arrow
        
        draw_quiver(x(root, 1), x(root, 2), x(leaf, 1), x(leaf, 2), msize, false);
        
    else % draw arrow at special points
        
        if sub_branch_num == 0 % reach an end of tree, big backgound
            draw_quiver(x(root, 1), x(root, 2), x(leaf, 1), x(leaf, 2), ...
                msize * 2, true, 'r');
        else % a branching point, medium background
            draw_quiver(x(root, 1), x(root, 2), x(leaf, 1), x(leaf, 2), ...
                msize * 2, true, 'y');
        end
        
        draw_quiver(x(root, 1), x(root, 2), x(leaf, 1), x(leaf, 2), msize, true);
        
    end % end if
    
end % end for 

end % end function