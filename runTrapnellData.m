%  runTrapnellData: Analysis Trapnell data with topographer
%
clear; clc; close all;
addpath(genpath(pwd))


% Step Controls; you can skip some steps by setting them to FALSE
DRAW_TRAJECTORY = 1 ;  % thow the trajectory of topographer
PLOT_TRANSITION = 1;  % transition probability & energy landscape
PLOT_PROFILE    = 1 ;  % plot gene profiles
DO_NETWORK      = 1;  % infer gene regulation network
DO_MODEL_INFER  = 1;  % infer expression parameters


%% data reading and filtering
dataName = 'trapnellRaw';
mFile = matfile([dataName '.mat']);
X_allrun = mFile.X;
X_raw = X_allrun(mFile.usefulCellIdx, :);
geneID_raw = mFile.geneID;

fprintf('Data from: %s \n', dataName);
fprintf('File descr: %s \n', mFile.description);
fprintf('Data scale (raw): N = %i, dim = %i \n', size(X_raw));
fprintf('Data sample: \n'); disp(X_raw(1:5, 1:5));

% Run a general filtering procedure
[X_prep, gene_ind] = exprPreprocessing(X_raw);
geneID_prep = geneID_raw(gene_ind);

% Check and exclude outliers (by clustering)
Z = linkage(X_prep, 'ward');
% dendrogram(Z); pause
clust_max = 20;
T = cluster(Z, 'maxclust', clust_max);
[pop_size, popCluster] = max(arrayfun(@(i) sum(T == i), 1:clust_max));
X_rm_outliers = X_prep(T == popCluster, :);
fprintf('removed %i/%i outliers.\n', numel(T) - pop_size, numel(T));

% Mamually remove some special cases
X_rm_outliers(234, :) = [];
X_rm_outliers(203, :) = [];
[~, X2d] = pca(X_rm_outliers, 'NumComponents', 2);
scatter(X2d(:, 1), X2d(:, 2), 30,'filled');
for i = 1:size(X2d, 1)
    text(X2d(i, 1), X2d(i, 2), num2str(i),'Color',[0.85 0.33 0.1]);
end

% Just use X and geneID for the rest of the script
X = X_rm_outliers;
geneID = geneID_prep;

% clear the confusing variables
clear mFile X_allrun X_raw geneID_raw X_prep geneID_prep X_rm_outliers X2d i
clear gene_ind cutoff clust_max Z T pop_size popCluster

fprintf('Data scale (filtered): N = %i, dim = %i \n', size(X));


%% (prepare for Topographer) (also prepare for gene-specific things)

% get some genes with high load in pca
coeff = pca(X); % Find some genes with high load in pca
[~, loadOrder] = sort(abs(coeff(:, 1)), 'descend');
highload_genes = loadOrder(1:4);
diff_expr_genes = DEG_find(X, 0);

% look for some marker genes
markersID = { ...
    'MYOG', 'MYF5' ... % TFs about myoblast differentiation
    'MYH2', ...        % marker for myotube matruation
    ...
    'CDK1', ...        % proliferation (early stage)
    'MEF2C', 'MYH2', ... % marker of (early / late) differentiation
    ...
    'PDGFRA', 'SPHK1', ... % not from myoblast lineage (contamination)
    };
markers_gene = zeros(1, numel(markersID));
for i = 1:numel(markersID)
    fprintf('Checking %s: ', markersID{i});
    fprintf('%i ', find(strcmpi(geneID, markersID{i})));
    fprintf('\n');
    markers_gene(i) = find(strcmpi(geneID, markersID{i}));
end

% what are the initial cells:
early_marker = find(strcmpi(geneID, 'CDK1'));
[~, early_order] = sort(X(:, early_marker), 'descend');
initial_cells = early_order(1:50);


%% topographer computation
opts.initialCells   = initial_cells;
opts.rQuantile      = 1;
opts.nTree          = 5*5;
opts.refineTimes    = 0;
opts.kManifold      = 15; % suggest:
opts.debug_all_traj = false;
opts.debug_refined  = false;
opts.debug_trimming = false;
opts.debug_print    = false;
opts.debug_get_traj = false;

S = topographer(X, opts);

% PCA for drawing
x_reduced = reduce2(X, 'pca');

% warning
if numel(S.tail_cells) == 1
    warning('Did not identify branching this time');
    disp('realy continue?');
    pause
end


%% Draw the final trajectory on 2D space
if DRAW_TRAJECTORY,
    figure(1);
    plotTrajectory(S, x_reduced);
end  % end if (step control)


%% Transition probability, and Energy landscape

if PLOT_TRANSITION,
    % calculate transition probability, by knn and pseudotime
    chi = 50; % time scaling factor
    verbose = false;
    [p_alpha, P, W] = directedTransition(S.pseudotime, S.knn_adj, chi, verbose);
    
    % group cell by state (using density cluster)
    LTree = linkedTree(S.trajectory_adj, S.progen_cell, S.t_assign); % tranform the trajectory into a nested tree structure
    groups_byState = groupByState(S.dist_mat, S.rho, LTree);
    
    % mannually (and interactively) pick and further divide the grouping
    figure(2);
    groups_byState = manDivide(groups_byState, S.pseudotime, x_reduced);
    
    % plot the transition
    figure(2);
    group_transition = plotTransition(P, groups_byState, x_reduced, S.rho);
    save_plot( 'trapnell_transition' , [], {'tif'});
    
    % define energy, by steady state probability
    energy = -log(p_alpha);
    energy = (energy - max(energy)) / range(energy); % normalize
    
    % plot energy landscape
    figure(3);
    gridNum = 50;
    smoothW = 3;  % smoothing window size
    moreF = true; % plot more information (cell location, trajectory)
    plot2DLandscape(x_reduced, S.trajectory_adj, energy, gridNum, smoothW, moreF);
    save_plot( 'trapnell_landscape' , [], {'tif'});
end % end if (step control)


%% Plot some gene profiles
if PLOT_PROFILE,
    figure(4);
    %plot_genes = highload_genes(1:3);
    plot_genes = markers_gene();
    plotGeneProfiles(S.pseudotime, X(:, plot_genes), geneID(plot_genes), LTree);
    save_plot('trapnell_Gene_Profiles', gcf, {'tif'});
end % end if (step control)


%% Peudotime-group-based gene analysis (GRN & inferModel)
% Group the cells by pseudotime, on a tree structure
time_nums = 5; % resulting in n groups
[groups_byTime, groupAdj] = groupByTime(LTree, ...
    S.pseudotime, time_nums, 'fuzz');

exprArray = cellfun(@(g) {X(g, :)} , groups_byTime'); % take group expression mat
% show group sizes
numGroup = numel(groups_byTime);
groups_size = arrayfun(@(i) length(groups_byTime{i}), 1:numGroup);
fprintf('time-group sizes: '); disp(groups_size);

% show the time grouping
figure(5);
debug_grouping(x_reduced, groups_byTime, groupAdj);

if DO_NETWORK
    
    collective = markers_gene(1:6);
    myoblast_marker = markers_gene(1:6);
    for baitGene = myoblast_marker(:)' % takes several minutes
        
        % find neighbour genes
        neighb_size = 5; % the 'k' parameter for kNN graph of genes
        min_conn = 2;
        mode = 'corr'; % genie mode is extremely slow !!!
        [neighbourGenes, ~] = networkNeighbor(X, baitGene, neighb_size, ...
            min_conn, mode);
        
        collective = union(collective, neighbourGenes);
    end
    
    % inferring regulations
    source = collective;
    target = source;
    verbose = true;
    GRNs = grnInference(exprArray, source, target, verbose);
    
    figure(6); clf;
    plotNetTree_simple(GRNs, groupAdj, geneID(source), 0.1);
    drawnow
    
    save_plot(sprintf('trapnell_GRN_of_markers_mode-%s', mode), [], {'fig', 'png'});
       
end % end if (step control)


if DO_MODEL_INFER
    
    infer_target = diff_expr_genes;
    warning('off')
    burst_param = parameterInference(exprArray, infer_target);
    warning('on')
    
end % end if (step control)
