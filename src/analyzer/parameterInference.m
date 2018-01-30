function plotNetTree_simple(grns, grn_net, node_id, weightTH)
% net net is a adj matrix representing relation among different GRNs

numNode = size(grns, 2);
numNet = size(grns, 1);

assert(size(grns, 2) == size(grns, 3));

if ~exist('node_id', 'var'), 
    node_id = arrayfun(@(i) {int2str(i)}, 1:numNode);
end
if ~exist('weightTH', 'var'), weightTH = 0; end

% circular layout (genes)
t = linspace(0, 2*pi, numNode + 1);
t = t(1:end-1);
x1 = cos(t); 
x2 = sin(t);
xnode = [x1', x2'];

% show each GRN
grnGrid = gridLayout(grn_net); % a naive grid-layout of tree nodes
[I, J, netIdx] = find(grnGrid);
[m, n] = size(grnGrid);
m = max(2, m); % handling the special case of no branching 
for iter = 1:numNet
    subplot(m, n, n * (I(iter)-1) + J(iter) );
    hold on
    
    gplot(squeeze(grns(netIdx(iter), :, :)) > weightTH, xnode);
    xlim([-1.1, 2]);
    ylim([-1.5, 1.5]);
    
    % maker the genes
    xtext = xnode;
    f_size = 8;
    for inode = 1:numNode
        if inode == 1
            text(xtext(inode, 1), xtext(inode, 2), node_id{inode}, ...
                'FontSize', f_size, 'FontWeight', 'bold', 'Color', 'r');
        else
            text(xtext(inode, 1), xtext(inode, 2), node_id{inode}, ...
                'FontSize', f_size);
        end
    end
    
    title( [ 'Stage ' int2str(netIdx(iter)) ] );
    axis off equal tight
    hold off
end

% show the relation between the tree and grid layout
subplot(m, n, n*(m-1)+1); % left-bottom position; usually not occupied
[~, stageOrder] = sort(netIdx);
coordi = nan(numel(I), 2);
coordi(:, 1) = J(stageOrder);
coordi(:, 2) = I(stageOrder);
gplot(grn_net, coordi, '-sr');
xlim([0.5, n + 0.5]); ylim([0.5, m+0.5]);
axis ij tight off

% step mark
hold on
for i = 1:numNet
    text(coordi(i, 1), coordi(i, 2), int2str(i), ...
        'FontSize', 12, 'FontWeight', 'bold');
end
hold off

drawnow

end

function mat = gridLayout(adj, root)
if ~exist('root', 'var'), root = find(~any(adj, 1)); end

leaf = find(adj(root, :));

if isempty(leaf),
    mat = root;
    return
end

mat = [];
for l = leaf(:)'
    mat = gridConcateVertical(mat, gridLayout(adj, l));
end
mat = gridConcateVertical(root, mat')'; % a trick to concate horizontally
end

function ret = gridConcateVertical(top, bot)
[m1, n1] = size(top);
[m2, n2] = size(bot);

ret = zeros(m1 + m2, max(n1, n2));

ret(1:m1, 1:n1) = top;
ret(m1+(1:m2), 1:n2) = bot;
end
