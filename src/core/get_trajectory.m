function S = get_trajectory(data, opts)
% finds the branching/tree-like trajectory multiplie times.
% The result is stored in structure S, along with some other usefull results.
% (see the end of this function for all returned results.)

% important global variables; working in elongate()
global T ridgecells

% The Algorithm Design
STRICT_RHO  = false; % true:  use the rho defined in Rodriguez' paper;
% false: use a gaussian kernal for rho.
LOCAL_WIDTH = false; % set a positvie value k: width of kernal <- dist to kth nn

% some alias
graph_k      = opts.kManifold;
traj_num     = opts.nTree;
d_c_quantile = opts.rQuantile;
db           = opts.debug_get_traj;

% some preparing calculations:
%   dist_mat is eulidean distance cumputed on the knn graph.
%   d_c is cutoff distance for local density; see Rodriguez's paper.
%   max_depth is the maximum depth of trajectory; set to prevent dead loop.
%   root_cell is the index of cell used as staring point for trajectories.
pre_tic = tic; % timing for debug printing
cell_num = size(data, 1);
[dist_mat, knn_adj] = graph_distance(data, graph_k, db);
dist_toc = toc(pre_tic);
d_c = quantile(nonzeros(triu(dist_mat)), d_c_quantile/100);
rho = local_density(dist_mat, d_c, STRICT_RHO, LOCAL_WIDTH);
max_depth = round(cell_num/10);
if isempty(opts.initialCells)
    [~, root_cell] = max(rho);
else
    root_cell = datasample(opts.initialCells,1);
end

if db,
    fprintf('\t preproc time: %.2f sec. ',  toc(pre_tic));
    fprintf('(%.2f sec for distance matrix) \n', dist_toc);
end

% find m trajectory tree, with randomized stepLen every step
tree_array = cell(1, traj_num);
for this = 1:traj_num
    traj_tic = tic;
    
    % NOTE: these global variables are initialized every run.
    % T is the adjacent matrix representing the found branching trajectory.
    % ridgecells will record all cells that are involved in the trajectory
    T = sparse(cell_num, cell_num);
    ridgecells = [];
    
    % THE recursive step
    depth = 1;
    pdt = dist_mat(root_cell, :); % an experimental try
    elongate(root_cell, [], depth, dist_mat, pdt, knn_adj, rho, d_c, max_depth, db);
    
    % store this trajactory tree
    tree_array{this} = T;
    
    if db,
        fprintf('\t %2d-th traj time: %.2f sec\n', this, toc(traj_tic));
    end
end

% record some results in S
S.dist_mat = dist_mat;       % distance matrix
S.rho = rho;                 % local density
S.d_c = d_c;                 % cutting distance; used as kernal width
S.knn_adj = knn_adj;         % adjacent matrix for the symmetric knn graph
S.root_cell = root_cell;     % cell index with the global density maximum
S.tree_array = tree_array;   % trajectory tree matrices (in a cellarray)

% clear the global working variables
clear global T ridgecells

end


function [D, knn_sym] = graph_distance(X, k, db)
% compute the pairwise distance between cells, using geodesic distance on a
% k-nearest-neighbour graph, witch is locally euclidean.

n = size(X, 1);

% determinate the k neighbors
[knnIdx, knnDist] = knnsearch(X, X, 'K', k + 1, 'Distance', 'euclidean');
knnIdx(:, 1) = []; % exclude self-link
knnDist(:, 1) = [];

% create the adjacent matrix
from = reshape(repmat(1:n, k, 1), 1, []);
to = reshape(knnIdx', 1, []);
w = reshape(knnDist', 1, []);
knn_weighted = sparse(from, to, w, n, n); % edges are distance-weighted
knn_sym = max(knn_weighted, knn_weighted'); % turn into undirected graph

D = graphallshortestpaths(tril(knn_sym), 'Directed', false);

if db && any(isinf(D))
    fprintf('warning: knn-graph not connected.');
    
    % incresase k to make a connected knn graph ?
    % fprintf('set k: %i -> %i\n', k, k+1);
    % k = k + 1;
    
    num_of_outliers = min(sum(isinf(D)));
    fprintf('Warning(get_trajectory/graph_distance): k manifold is not connected. ');
    fprintf('outlied: %d \n', num_of_outliers);
end

D = (D + D') / 2; % D is slightly unsymmetric due to computation error

end


%%%
% enlongate() is a very important (and long) function
function elongate(this, prev, depth, dist_mat, pdt, knn, rho, d_c, max_depth, db)
% A recursive function, modifying global variables 'T' and 'ridgecells'.

% they are initialized in the main function.
global T ridgecells

% The Algorithm Design
DELTA_ONLY = false;  % true: return K points with top delta value
LOCAL_DELTA = false;  % true: re-calculate density on the ring. SLOW
K = 3;

% some parameters
step_factor = 2.5;  % step length, relative to local sigma
ring_thick = 1/5;   % 1/5 of the 2.5-sig step -> ring = 2-sig ~ 3-sig
variation = 0.1;    % introduce how much random variation in step len


%% Self check to prevent collision

if depth > max_depth
    if db, % only print wranning for debug
        disp('Warning(get_trajectory/elongate): exceed maxDepth. forced termination.');
    end
    return;
end

if depth > 1
    % debug checking
    assert(ismember(prev, ridgecells));
    
    % not going back
    if ismember(this, ridgecells), return; end
    if pdt(this) < pdt(prev), return; end
    
    % not too close to other ridge cells
    is_neighbor_of_others = any(ismember(find(knn(this, :)), ridgecells));
    if is_neighbor_of_others, return; end
    
    oldercels = ridgecells(ridgecells ~= prev);
    shared_neighbor_with_oldercells = any( knn(this, :) & any(knn(oldercels, :), 1) );
    if shared_neighbor_with_oldercells, return; end
end

% pass the check, so put this into the trajectory
ridgecells = [ridgecells this];
T(prev, this) = true; % add 'this' cell into the collection of ridges


%% Find peaks, and filter out bad peaks

% a shortcut function
sigma = @(i) ridge_width(i, dist_mat, rho, d_c);
step_len = sigma(this) * step_factor * (1 + variation * randn); % randomized

% find peak by gamma
while true % loop for sweeping search, only at the first step
    try
        peaks = find_peaks(this, step_len, dist_mat, rho, ring_thick, K, DELTA_ONLY, LOCAL_DELTA);
        break
    catch ME
        if strcmp(ME.identifier, 'topographer:noPeaks')
            if depth == 1
                dist2ring = abs(dist_mat(this, :) - step_len);
                [~, targetCell] = min(dist2ring);
                
                % change the stepLen if no peaks found
                if step_len < dist_mat(this, targetCell), act = 'Expand';
                else act = 'Shrink'; end
                
                if db,  % only print warnning for debug
                    fprintf('Warning(get_trajectory/elongate): %s stepLen for sweepingSearch.', act);
                end
                
                step_len = dist_mat(this, targetCell);
                continue % try again
            else
                % reached an terminal. just stop the search
                return
            end
        else
            rethrow(ME);
        end
    end % end try
end % end loop


%% Resursive part, going deeper
% NOTE: nothing happens if 'peaks' is [] !

peaks = sort_by(peaks, rho(peaks), 'descend');
for p = peaks(:)'
    depth = depth + 1;
    elongate(p, this, depth, dist_mat, pdt, knn, rho, d_c, max_depth, db);
end


end


%%%
% Some helper functions for enlongate()

function r = ridge_width(x, dist, rho, d_c)
% estimating the width of moutain ridge, assuming gaussian shape.

global DEBUG_SIGMA % for debug

if length(x) > 1 % a little trick to handel multiple inputs
    r = arrayfun( @(xi) ridge_width(xi, dist, rho, d_c) , x);
else
    % Find the closest cell x such that rho(x) < sigmaRho
    [~, distOrder] = sort(dist(x, :), 'ascend');
    sigma_rho = rho(x) * 0.606531; % 0.6065 = 1/sqrt(e) = norm(sigma)/norm(0)
    sigma_cell = distOrder(find(rho(distOrder) <= sigma_rho, 1, 'first'));
    if isempty(sigma_cell)
        r = d_c;
    else
        sigma_dist = dist(x, sigma_cell);
        r = min(sigma_dist, 2 * d_c); % prevent it from geting too large
        
        DEBUG_SIGMA(x, sigma_cell) = sigma_dist; % record for debug
    end
end

end


function peaks = find_peaks(x, radius, D, rho, width, K, deltaOnly, localData)
% findPeaks finds gamma peaks on the ring using Rodriguez's method.
% x: center of the ring. an index.

ringPoints = find_points_on_ring(D(x, :), radius, width);
%assert(length(ringPoint) >= K, 'valle:noPeaks', 'Not enough points');
assert(length(ringPoints) >= 1, 'topographer:noPeaks', 'Not enough points');

if localData
    ringD = D(ringPoints, ringPoints);
    ringRho = local_density(ringD, radius/2);
else
    ringD = D(ringPoints, ringPoints);
    ringRho = rho(ringPoints);
end

peakIdx = gamma_peak(ringD, ringRho, min(K, length(ringPoints)), deltaOnly);
peaks = ringPoints(peakIdx);

end


function ringIndex = find_points_on_ring(dist2center, radius, width)
innerR = radius * (1 - width);
outerR = radius * (1 + width);
ringIndex = find( dist2center>=innerR & dist2center<=outerR );
end


function peaks = gamma_peak(dist, rho, k, useDeltaOnly)
%GAMAPEAK call density_peaks() to find density peaks. The results are
% ordered in 'gamma'.
%
% NOTE: if useDeltaOnly is true, gammaPeak would return exactly k clusters.

[g_peaks, delta, gamma] = density_peaks(dist, rho, k);

if useDeltaOnly
    
    % use delta instead of gamma
    [~, delta_order] = sort(delta, 'descend');
    peaks = delta_order(1:k);
    
else
    
    % order the results in gamma
    peaks = sort_by(g_peaks, gamma(g_peaks), 'descend');
    
end

end