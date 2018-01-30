function S = get_refined_tree(S, ~, opts)
%REFINETREE do clustering on the treeArray, and pick the most representive
% trajectory tree as S.T. This S.T is also refined by grahpInsert and
% graphExtend.
%
% @NOTE: since pdt is used in getTrajectory(), trimming is useless now. 

db = opts.debug_trimming;
trimming = false;
refining = true;

% alias
knn = S.knn_adj;
dist = S.dist_mat;
root = S.root_cell;
tree_array = S.tree_array;
pdt = dist(root, :);

% trim each found trajectory
if trimming,
    for i = 1:opts.nTree
        T = tree_array{i};
        T = graphInsert(T, dist, S.rho, opts.insertWidth);
        
        if db, fprintf('- %dth: \n', i); end
        T = trimBranchRecur(T, knn, dist, pdt, root, db);
        
        tree_array{i} = T; % update
    end
end

% select the trajectory with best coverage
if opts.nTree > 1,
    %[S, select] = treeClustering(S);
    
    % compute coverage
    getTreeNodes = @(t) [find(any(t(:, :), 1)) root];
    getCovered = @(ids) union(find(any(knn(ids, :), 1)), ids);
    getCoverage = @(t) numel(getCovered( getTreeNodes(t) ));
    treeCoverage = cellfun(getCoverage, tree_array);
    [~, select] = max(treeCoverage);
    
else
    select = 1;
end

% 5. do the "insert" and "extend" refinement
M = tree_array{select};
if refining,
    for i = 1:opts.refineTimes,
        M = graphInsert(M, dist, S.rho, opts.insertWidth);
        %M = graphExtend(M, S.D, S.rho, opts.extDegree);
    end
end % end if

S.trimmed_tree_array = tree_array;  
S.selected_tree_idx = select;  % the index of finally selected tree
S.refined_tree = M;

end


% Tree clustering (OLD, not used now)
function [S, select] = treeClustering(S)
% just a place to put the old code for tree clustering

% 1. calcualte pair-wise distance of Trees
treeD = zeros(opts.nTree);
for i = 2:opts.nTree
    for j = 1:(i-1)
        d = treeDist(S.treeArray{i}, S.treeArray{j}, S.D);
        treeD(i, j) = d;
        treeD(j, i) = d;
    end
end
% 2. do hierachial clustering
if opts.nTreeCluster > 1
    z = linkage(squareform(treeD), 'average');
    c = cluster(z, 'maxclust', opts.nTreeCluster);
else
    c = ones(1, opts.nTree);
    z = [];
end
S.treeLinkage = z;      % the dendrogram of the hierachial clustering
S.treeCluster = c;      % the tree clustering result

% 3. find representative "centroid" for each cluster,
nCluster = max(c);
centroidTrees = zeros(nCluster, 1);
centroidCost = zeros(nCluster, 1);
for icluster = 1:nCluster
    % find the "centroid" in this cluster
    clusterMember = find(c == icluster);
    [~, centroidIndex] = min( sum(treeD(clusterMember, clusterMember).^ 2));
    icentroid = clusterMember(centroidIndex);
    
    % calculate min distance for each cell to the tree
    [I, J] = find( S.treeArray{icentroid} );
    treeNode = unique([I J]);
    [minDist, ~] = min(S.D(:, treeNode), [], 2);
    
    centroidTrees(icluster) = icentroid;
    centroidCost(icluster) = sum(minDist .^ 2) ^ 0.5;
end

% 4. select the centroids with best coverage,
%    and a tiny consideration on cluster scale ('memberCounts')
memberCounts = arrayfun(@(i) sum(c == i), 1:nCluster);
[~, selectIndex] = min(centroidCost .^ 2 ./ memberCounts');
select = centroidTrees(selectIndex);
end
function d = treeDist(T1, T2, D)
%TREEDIST defines the distance between two trees

% make sure the matrix are not symmetric
assert(~issymmetric(T1) && ~issymmetric(T2));

% get the edges in the graphs
E1 = getEdges(T1);
E2 = getEdges(T2);

% calculare pair-wise edge distance
edgeD = inf(size(E1, 1), size(E2, 1)); % initialization
for i = 1:size(E1, 1)
    for j = 1:size(E2, 1)
        edgeD(i, j) = edgeDist(E1(i, :), E2(j, :), D);
    end
end
assert(~any(isnan(nnz(edgeD)))); % check for bug

% set d to be the sum of edge-to-graph distance
d1 = sum(min(edgeD, [], 2) .^ 2 );
d2 = sum(min(edgeD, [], 1) .^ 2 );
d = max(d1, d2);
end
function e = getEdges(T)
[I, J] = find(T);
e = [I J];
end
function d = edgeDist(a, b, D)
% define distance bwtween two edge.
d1 = D(a(1), b(1)) + D(a(2), b(2));
d2 = D(a(1), b(2)) + D(a(2), b(1));
d = min(d1, d2);
end


% Tree refining

function M = graphExtend(M, D, rho, degree)
% GRAPHEXTEND extend the ending segment of the graph
% note: M is assumed to be an underected graph matrix

assert(issymmetric(M));
endNodes = find(sum(M, 1) == 1);
for i = 1:length(endNodes)
    p = endNodes(i);
    prev = find(M(p, :)); assert(length(prev) == 1);
    
    % find the points ahead, in a range defined by 'degree'
    compleAng = (180 - degree) / 180 * pi;
    inRange = @(x) acosByDist(D(p,prev), D(x,p), D(x,prev)) > compleAng;
    isInRange = arrayfun(@(x) inRange(x), 1:length(rho));
    candidate = find(isInRange);
    
    % prevent extend to points too far
    step = D(p, prev);
    candidate = candidate(D(p, candidate) > step & D(p, candidate) < 2*step);
    % make sure that no loop wiil be  produced
    % this can be canceled in kManifold
    
    if isempty(candidate), continue; end
    
    % select the candidate cell with highest local density
    [~, tmp] = max(rho(candidate));
    ext = candidate(tmp);
    
    % update M
    M(p, ext) = true; M(ext, p) = true;
end

end

function M = graphInsert(M, D, rho, width)
% REFINEFRAPH insert immediate point at the edge, to refine the graph.
% note: M is assumed to be an directed, acyclic simple graph matrix; if M is 
% symmetric, we use only triu(M).

sym_flag = false;

if issymmetric(M)
    sym_flag = true;
    disp('Warning (graphInsert): symmtric graph is inputed.');
else
    assert(~any(reshape(M & M', 1, []))); % make sure links to be simple
end

% insert intermedian points
[I, J, W] = find(M);
ridgeCells = union(I', J');
checkList = sortBy(1:length(W), rho(I) + rho(J), 'descend');
for k = checkList(:)'
    i = I(k); j = J(k);
    w = W(k);
    
    r = D(i, j) * (0.5 + width);
    candidate = find( (D(:, i) <= r) & (D(:, j) <= r) );
    candidate((candidate == i) | (candidate == j)) = []; % prevent error
    candidate(ismember(candidate, ridgeCells)) = []; % prevent loop
    if isempty(candidate), continue; end
    
    % select the cell with highest density
    [~, ind] = max(rho(candidate));
    inter = candidate(ind);
    
    % update M
    M(i, j) = 0; % break the old edge
    M(i, inter) = w;
    M(inter, j) = w;
    ridgeCells = [ridgeCells inter];
end

if sym_flag
    M = M + M'; % make output symmetric as the input
end

end

function T = trimBranchRecur(T, knn, dist, pdt, root, db)
leaf = find(T(root, :));
assert(~isempty(leaf));

% 1. trim back-heading branches
% check back-headding if branches are simple (or just the backbone part)
% filter based on directed-direction (matuarity)
minAdvance = 0;
for l = leaf(:)'
    prev = root;
    this = l;
    while numel(this) == 1 % make sure it is on a simple branch
        if pdt(this) < pdt(prev) + minAdvance,
            % find the end of 'back-heading'
            jump = find(T(this, :));
            while (numel(jump) == 1) && ...
                    (pdt(jump) < pdt(prev) + minAdvance)
                jump = find(T(jump, :));
            end
            if numel(jump) ~= 1
                % end of bakcbone -> trim all the rest
                T = cleanDown(T, this);
                if db, disp('backheading: trim a branch'); end
                break; % end searching
            else
                % skip the intermediate part by jumping
                T = cleanDown(T, this, jump);
                T(prev, jump) = true;
                if db, disp('backheading: skip a segment'); end
                prev = jump; % jumping
                this = find(T(prev, :)); % update, since 'prev' is changed
            end
        else
            % just step ahead
            prev = this;
            this = find(T(prev, :)); % update, since 'prev' is changed
        end
    end
end

% 2. merge the branches, if there are two
getCover = @(ids) union(find(any(knn(ids, :), 1)), ids); % useful
getCoverN = @(ids) numel(getCover(ids));
get2ndCover = @(ids) getCover(getCover(ids)); % useful
if numel(leaf) == 1,
    % only one branch. no nead to merge. just end.
    return
end
% sort the branches by subtree cover size. we check the smaller branch first.
subtreeCoverSizes = arrayfun(@(l) getCoverN(getTreecell(l, T, 'all')), leaf);
leaf = sortBy(leaf, subtreeCoverSizes, 'ascend');
% pair-wise merge. 
for i = 1:(numel(leaf)-1)
    l = leaf(i);
    iBackbone = getTreecell(l, T, 'backbone');
    
    % check if the i-th branch is adhered to a bigger branch
    for j = (i+1):length(leaf)
        jCells = getTreecell(leaf(j), T, 'all');
        detachIdx = find(~ismember(iBackbone, get2ndCover(jCells)), 1, 'first');
        if isempty(detachIdx)
            % no detachment -> trim the entire subbranch
            T = cleanDown(T, l);
            if db, disp('branching: entire trimming'); end
        elseif detachIdx > 1
            % trim part of 'iBackbone' that is adhered to jCells
            holder = iBackbone(detachIdx);
            T = cleanDown(T, l, holder); % trimming happens here
            % reconnect
            lastAdhered = iBackbone(detachIdx-1);
            [~, closestIdx] = min(dist(lastAdhered, jCells));
            T(jCells(closestIdx), holder) = true;
            if db, disp('branching: re-connect'); end
        else
            % no need to trim
            continue
        end
        
        % trimming happens, so it is done for the i-th branch
        break;
    end
end

% Finally, take care of subtree in all branches (if any) recursively.
leaf = find(T(root, :)); % update, bacause it may has been changed
for l = leaf(:)'
    subroot = l;
    next = find(T(subroot, :));
    while numel(next) == 1;
        subroot = next;
        next = find(T(subroot, :));
    end
    if numel(next) > 1
        T = trimBranchRecur(T, knn, dist, pdt, subroot, db);
    end
end
end


% little helpers

function branchCells = getTreecell(root, T, mode)
% result include 'root'

leaves = find(T(root, :));
if isempty(leaves),
    branchCells = root;
    return
end

if length(leaves) == 1
    branchCells = [root getTreecell(leaves, T, mode)];
else % have multiple
    switch mode
        case 'backbone'
            branchCells = root;
        case 'all'
            branchCells = root;
            for l = leaves(:)',
                branchCells = [branchCells getTreecell(l, T, mode)];
            end
    end
end

end

function T = cleanDown(T, start, stop)
% Remove nodes in the graph T. Include 'start', possiblely not include 'stop'.
if ~exist('stop', 'var'), stop = -1; end

if start == stop, return; end

T(:, start) = false;
for l = find(T(start, :));
    T = cleanDown(T, l, stop);
end
T(start, :) = false;
end

function result = sortBy(x, orderdata, mode)
[~, order] = sort(orderdata, mode);
result = x(order);
end

function ang = acosByDist(a, b, c)
ang = acos( (a^2 + b^2 - c^2) / (2 * a * b) );
ang = real(ang); % some time be complex, for non-euclidean input
end