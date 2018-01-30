function S = get_branches(S, data, opts)

% alias and preprocessing 
undirected_traj = min(S.refined_tree + S.refined_tree', 1);
initial_cells = opts.initialCells;


%% 1. select a progenitor cells on the trajectory, using initial cells

node_degree = sum(full(undirected_traj));
end_nodes = find(node_degree == 1);   % used when no initial cells proveide 

traj_cells = find(node_degree > 0); % all cells on the trajectory
if isempty(initial_cells)
    disp('Warning (gerBranches): no initial cells provided; use random cells.');
    progen_cell = randsample(end_nodes, 1);
else
    assert(~any(isinf(S.dist_mat( initial_cells, S.root_cell )))); % debug test 
    [~, idx] = min(arrayfun(@(t) sum(S.dist_mat(t, initial_cells)), traj_cells));
    progen_cell = traj_cells(idx);
end


%% 2. reidrecte the graph using the selected progenitor cell (progen_cell)

% 1. generate paths by Breath-First Searching, through a built-in function
%    using this function would be inefficient; but convinient for now. 
[~, ~, pred] = graphshortestpath(undirected_traj, progen_cell, 'Method', 'BFS');
% generate a directed trajAdj, coresponsed to M
mark = pred > 0; % mark all leaf nodes; nan > 0 is also false
from = pred(mark);
to = find(mark);
n = size(undirected_traj, 1);
trajectory_mat = sparse(from, to, ones(length(from), 1), n, n);

% check: they should be identical, aside from direction 
assert(all(all(undirected_traj == min(trajectory_mat + trajectory_mat', 1))));

% 2. determine tail nodes
tail_cells = end_nodes(end_nodes ~= progen_cell);

% % 3. get branches from 'path'
% nBranches = length(tailing_cells);
% branchArray = cell(nBranches, 1); % initialization
% for i = 1:nBranches
%     branchArray{i} = path{ tailing_cells(i) };
% end
% S.nBranches = nBranches;
% S.branchArray = branchArray;
% S.predecessor = pred;


%% 3. Assign all other cells to one of the trajectory cells

% for those cell not connected in the k manifold, link them to a closest trajector cell
disconnected_cells = find(isinf(S.dist_mat(progen_cell, :)));
[closest_idx, d] = knnsearch(data(traj_cells, :), data(disconnected_cells, :), 'K', 1);
closest_traj_cell = traj_cells(closest_idx);
for i = 1:length(disconnected_cells)
    S.dist_mat(disconnected_cells(i), closest_traj_cell(i)) = d(i);
    S.dist_mat(closest_traj_cell(i), disconnected_cells(i)) = d(i);
end

% assign all cells to a closest trajectory cell
[~, assign_idx] = min(S.dist_mat(:, traj_cells), [], 2);
traj_assign = traj_cells(assign_idx);


%% return

S.progen_cell = progen_cell;
S.trajectory_adj = trajectory_mat;
S.trajectory_pred_link = pred;
S.trajectory_cells = traj_cells;
S.tail_cells = tail_cells;
S.t_assign = traj_assign;


end