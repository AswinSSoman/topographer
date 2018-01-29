% runTestData: Test the topographer algorithm on synthetic data.
%
clear; clc; close all;
addpath(genpath(pwd))


%% Select a Data Types
option                  = {'S', 'Y', 'DoubleBranch', 'Spiral', 'Wadd', 'Triple'};
dataType                = option{ 5 };
nSamples                = 1000;
[tdata, theFirstPoints] = testData(dataType,nSamples);
X                       = tdata;

% draw the data set for debug 
% figure
% scatter(X(:, 1), X(:, 2), 30);


%% Topographer computation   
opts.initialCells   = theFirstPoints;
opts.rQuantile      = 2;
opts.nTree          = 9; 
opts.refineTimes    = 0;
opts.kManifold      = 30; 
opts.debug_all_traj = false;
opts.debug_refined  = false;
opts.debug_trimming = false;
opts.debug_print    = true;
opts.debug_get_traj = false;

S = topographer(X, opts);  % see the content of S 


%% Draw the final trajectory on 2D space  
% PCA for drawing 
x_reduced = reduce2(X, 'pca');
figure(1);
plotTrajectory(S, x_reduced);


%% Transition probability, and Energy landscape   

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
figure(3);
group_transition = plotTransition(P, groups_byState, x_reduced, S.rho);
save_plot( 'test_transition' , [], 'tif');

% define energy, by steady state probability
energy = -log(p_alpha);
energy = (energy - max(energy)) / range(energy); % normalize

% plot energy landscape 
figure(4);
gridNum = 50; 
smoothW = 3;  % smoothing window size
moreF = true; % plot more information (cell location, trajectory) 
plot2DLandscape(x_reduced, S.trajectory_adj, energy, gridNum, smoothW, moreF);
save_plot( 'test_landscape' , [], 'tif');

