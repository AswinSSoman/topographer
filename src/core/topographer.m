function S = topographer(data, input_options)
%TOPOGRAPHER Carray out a seires of computation to find the trajectory in data
% return a struct containing all usefull results.
%
% @param data is a m-by-n matrix, each row is expression vector for a cell. 
% @param input_opionts is a struct specifying the parameters for topographer.
% input_opionts.initialCells is required for meaningful computation. 

global SIGMA  % DEBUG variable 

%% Algorithm defaults

opts.initialCells   = [];       % must be provided by user
opts.rQuantile      = 0.2;      % default in the gamma-cluster paper
opts.nTree          = 15;       % number of trajectory trees to find
opts.nTreeCluster   = 4;        % number of tree clusters
opts.extDegree      = 30;       % range to extend from end node of the graph
opts.insertWidth    = 0.15;     % the portion of an edge to insert an extra node
opts.refineTimes    = 1;        % times of repeating refinement
opts.kManifold      = 10;

opts.debug_all_traj = 1;
opts.debug_refined  = false;
opts.debug_trimming = false;
opts.debug_print    = false;
opts.debug_get_traj = false;

% Read the input options
intpu_params = fieldnames(input_options);
for i = 1:length(intpu_params)
    name = intpu_params{i};
    value = input_options.(name);
    if isfield(opts, name);
       opts.(name) = value;
    else
       fprintf('Warning(topographer): invalid option "%s"', name);
    end
end

% initialization
func_tic = tic;
n = size(data, 1);
SIGMA = sparse(n, n); % only used for debug purpose 


%% 1. Get the trajectory
% Function getTrajectory() finds possible trajectory tree multiple times, and 
% stores the result in S.treeArray. Each tree in S.treeArray is a
% adjacent matrix representing a simple acyclic directed graph, rooted at the
% cell with highest local density.

tic;
S = get_trajectory(data, opts);
if opts.debug_print, fprintf('traj time: %.2f sec\n', toc); end


if opts.debug_all_traj
    tic; figure(98); clf;
    debug_all_traj(S, data, opts);
    fprintf('dbat time: %.2f sec\n', toc);
end


%% 2. refine: clustering?, extend, insert
% most important step in this stage is filtering all found trees, and find out 
% the (or synthesize a) most representative one. The result is stored in 
% S.refined_tree. 

tic;
S = get_refined_tree(S, data, opts);
if opts.debug_print, fprintf('refn time: %.2f sec\n', toc); end


%% 3. map cells on branches
% some of the most complicated (but less important) computation happens here. 
% getBranches() sorts the results of previous steps and returns a more useful
% tructure for later computation. The key step is redirecting the
% S.refined_tree using the user-provided initial cell(s). The resulted tree is
% then a more meaningfull branching trajectory, beginning from a plausible 
% progenitor cell ( instead of a less meaningful, local-density-maximum cell ).

tic;
S = get_branches(S, data, opts);
if opts.debug_print, fprintf('brnc time: %.2f sec\n', toc); end


if opts.debug_refined
    tic; figure(99); clf;
    debug_selecting(S, data, 'pca');
    fprintf('dbpl time: %.2f sec\n', toc);
end


%% 4. calculate Pseudotime
tic;
S = get_pseudotime(S, data, opts);
if opts.debug_print, fprintf('psdt time: %.2f sec\n', toc); end

if opts.debug_print, 
    fprintf('Topographer done. time: %.2f sec\n', toc(func_tic)); 
end