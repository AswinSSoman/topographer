function fig = wdlst(X, source, geneToshow, geneNamesArray, fig)
  %WDLST is just a wapper for the wanderlust code, calling wanderlust.m and
  %    wanderlust_line_graph.m
  %
  
  % run the wanderlust algorithm and get the trajectorys
  data = X;
  distance = 'euclidean';
  k = 8; l = 11;
  nGraphs = 15;
  nLandmarks = 10;
  verbose = true;
  partial_order = []; % list of already ordered cells
  [graph_traj, lnn] = wanderlust( data, distance,...
    k, l, nGraphs, source, nLandmarks, verbose, partial_order);
  
  %   % show the lnn graph
  %   figure(fig); fig = fig +1; clf;
  %   if size(X, 2) > 2 % tsne
  %     pca_dims = min(size(X, 2), 50); % tSNE applys PCA as preprocessing
  %     perplexity = 30; % default is 30
  %     X2d = tsne(X, zeros(1, size(X, 1)), 2, pca_dims, perplexity);
  %   else
  %     X2d = X;
  %   end
  %   gplot(lnn, X2d,  '--*'); hold on
  %   scatter(X2d(source, 1), X2d(source, 2), 100, 'k'); % circle the source
  
  nGene = length(geneToshow);
  n = ceil(nGene^0.5);
  m = ceil(nGene/n);
  
  % plot the wanderlust line
  wandelust_reaulst = mean(graph_traj);
  num_pts = 15; % time step numbers
  window_width = 0.1;
  figure(fig); fig = fig +1; clf;
  for i = 1:nGene
    gene = geneToshow(i);
    name = geneNamesArray{i};
    
    data = X(:, gene);
    groupAssign = [];
    groupLabel = {''};
    
    subplot(m, n, i);
    [ pts, pt_indices, pt_mean_intensity, pt_std_intensity ] = ...
      wanderlust_line_graph( wandelust_reaulst, data, groupLabel, num_pts, ...
      window_width, groupAssign );
    title(name);
  end
end