function ydata = tsne(X, labels, no_dims, initial_dims, perplexity, verbose)
  %TSNE Performs symmetric t-SNE on dataset X
  %
  %   mappedX = tsne(X, labels, no_dims, initial_dims, perplexity)
  %   mappedX = tsne(X, labels, initial_solution, perplexity)
  %
  % The function performs symmetric t-SNE on the NxD dataset X to reduce its
  % dimensionality to no_dims dimensions (default = 2). The data is
  % preprocessed using PCA, reducing the dimensionality to initial_dims
  % dimensions (default = 30). Alternatively, an initial solution obtained
  % from an other dimensionality reduction technique may be specified in
  % initial_solution. The perplexity of the Gaussian kernel that is employed
  % can be specified through perplexity (default = 30). The labels of the
  % data are not used by t-SNE itself, however, they are used to color
  % intermediate plots. Please provide an empty labels matrix [] if you
  % don't want to plot results during the optimization.
  % The low-dimensional data representation is returned in mappedX.
  %
  %
  % (C) Laurens van der Maaten, 2010
  % University of California, San Diego
  
  % thie variable is added by syaoming
  if ~exist('verbose', 'var'), verbose = false; end
  
  if ~exist('labels', 'var')
    labels = [];
  end
  if ~exist('no_dims', 'var') || isempty(no_dims)
    no_dims = 2;
  end
  if ~exist('initial_dims', 'var') || isempty(initial_dims)
    initial_dims = min(50, size(X, 2));
  end
  if ~exist('perplexity', 'var') || isempty(perplexity)
    perplexity = 30;
  end
  
  % First check whether we already have an initial solution
  if numel(no_dims) > 1
    initial_solution = true;
    ydata = no_dims;
    no_dims = size(ydata, 2);
    perplexity = initial_dims;
  else
    initial_solution = false;
  end
  
  % Normalize input data
  X = X - min(X(:));
  X = X / max(X(:));
  X = bsxfun(@minus, X, mean(X, 1));
  
  % Perform preprocessing using PCA
  if ~initial_solution
    if verbose, disp('Preprocessing data using PCA...'); end
    if size(X, 2) < size(X, 1)
      C = X' * X;
    else
      C = (1 / size(X, 1)) * (X * X');
    end
    [M, lambda] = eig(C);
    [lambda, ind] = sort(diag(lambda), 'descend');
    M = M(:,ind(1:initial_dims));
    lambda = lambda(1:initial_dims);
    if ~(size(X, 2) < size(X, 1))
      M = bsxfun(@times, X' * M, (1 ./ sqrt(size(X, 1) .* lambda))');
    end
    X = bsxfun(@minus, X, mean(X, 1)) * M;
    clear M lambda ind
  end
  
  % Compute pairwise distance matrix
  sum_X = sum(X .^ 2, 2);
  D = bsxfun(@plus, sum_X, bsxfun(@plus, sum_X', -2 * (X * X')));
  
  % Compute joint probabilities
  P = d2p(D, perplexity, 1e-5, verbose); % compute affinities using fixed perplexity
  clear D
  
  % Run t-SNE
  if initial_solution
    ydata = tsne_p(P, labels, ydata, verbose);
  else
    ydata = tsne_p(P, labels, no_dims, verbose);
  end
end


function ydata = tsne_p(P, labels, no_dims, verbose)
  %TSNE_P Performs symmetric t-SNE on affinity matrix P
  %
  %   mappedX = tsne_p(P, labels, no_dims)
  %
  % The function performs symmetric t-SNE on pairwise similarity matrix P
  % to create a low-dimensional map of no_dims dimensions (default = 2).
  % The matrix P is assumed to be symmetric, sum up to 1, and have zeros
  % on the diagonal.
  % The labels of the data are not used by t-SNE itself, however, they
  % are used to color intermediate plots. Please provide an empty labels
  % matrix [] if you don't want to plot results during the optimization.
  % The low-dimensional data representation is returned in mappedX.
  %
  %
  % (C) Laurens van der Maaten, 2010
  % University of California, San Diego
  
  
  if ~exist('labels', 'var')
    labels = [];
  end
  if ~exist('no_dims', 'var') || isempty(no_dims)
    no_dims = 2;
  end
  
  % First check whether we already have an initial solution
  if numel(no_dims) > 1
    initial_solution = true;
    ydata = no_dims;
    no_dims = size(ydata, 2);
  else
    initial_solution = false;
  end
  
  % Initialize some variables
  n = size(P, 1);                                     % number of instances
  momentum = 0.5;                                     % initial momentum
  final_momentum = 0.8;                               % value to which momentum is changed
  mom_switch_iter = 250;                              % iteration at which momentum is changed
  stop_lying_iter = 100;                              % iteration at which lying about P-values is stopped
  % max_iter = 1000;                                    % maximum number of iterations
  max_iter = 300; % mod by syaoming
  epsilon = 500;                                      % initial learning rate
  min_gain = .01;                                     % minimum gain for delta-bar-delta
  
  % Make sure P-vals are set properly
  P(1:n + 1:end) = 0;                                 % set diagonal to zero
  P = 0.5 * (P + P');                                 % symmetrize P-values
  P = max(P ./ sum(P(:)), realmin);                   % make sure P-values sum to one
  const = sum(P(:) .* log(P(:)));                     % constant in KL divergence
  if ~initial_solution
    P = P * 4;                                      % lie about the P-vals to find better local minima
  end
  
  % Initialize the solution
  if ~initial_solution
    ydata = .0001 * randn(n, no_dims);
  end
  y_incs  = zeros(size(ydata));
  gains = ones(size(ydata));
  
  % Run the iterations
  for iter=1:max_iter
    
    % Compute joint probability that point i and j are neighbors
    sum_ydata = sum(ydata .^ 2, 2);
    num = 1 ./ (1 + bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata', -2 * (ydata * ydata')))); % Student-t distribution
    num(1:n+1:end) = 0;                                                 % set diagonal to zero
    Q = max(num ./ sum(num(:)), realmin);                               % normalize to get probabilities
    
    % Compute the gradients (faster implementation)
    L = (P - Q) .* num;
    y_grads = 4 * (diag(sum(L, 1)) - L) * ydata;
    
    % Update the solution
    gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...         % note that the y_grads are actually -y_grads
      + (gains * .8) .* (sign(y_grads) == sign(y_incs));
    gains(gains < min_gain) = min_gain;
    y_incs = momentum * y_incs - epsilon * (gains .* y_grads);
    ydata = ydata + y_incs;
    ydata = bsxfun(@minus, ydata, mean(ydata, 1));
    
    % Update the momentum if necessary
    if iter == mom_switch_iter
      momentum = final_momentum;
    end
    if iter == stop_lying_iter && ~initial_solution
      P = P ./ 4;
    end
    
    % Print out progress
    if verbose && ~rem(iter, 50) % mod by syaoming, used to be ~rem(iter, 10)
      cost = const - sum(P(:) .* log(Q(:)));
      disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]);
    end
    
    % Display scatter plot (maximally first three dimensions)
    if ~rem(iter, 10) && ~isempty(labels)
      if no_dims == 1
        scatter(ydata, ydata, 9, labels, 'filled');
      elseif no_dims == 2
        scatter(ydata(:,1), ydata(:,2), 9, labels, 'filled');
      else
        scatter3(ydata(:,1), ydata(:,2), ydata(:,3), 40, labels, 'filled');
      end
      axis tight
      axis off
      drawnow
    end
  end
end

function ydata = tsne_d(D, labels, no_dims, perplexity)
  %TSNE_D Performs symmetric t-SNE on the pairwise Euclidean distance matrix D
  %
  %   mappedX = tsne_d(D, labels, no_dims, perplexity)
  %   mappedX = tsne_d(D, labels, initial_solution, perplexity)
  %
  % The function performs symmetric t-SNE on the NxN pairwise Euclidean
  % distance matrix D to construct an embedding with no_dims dimensions
  % (default = 2). An initial solution obtained from an other dimensionality
  % reduction technique may be specified in initial_solution.
  % The perplexity of the Gaussian kernel that is employed can be specified
  % through perplexity (default = 30). The labels of the data are not used
  % by t-SNE itself, however, they are used to color intermediate plots.
  % Please provide an empty labels matrix [] if you don't want to plot
  % results during the optimization.
  % The data embedding is returned in mappedX.
  %
  %
  % (C) Laurens van der Maaten, 2010
  % University of California, San Diego
  
  
  if ~exist('labels', 'var')
    labels = [];
  end
  if ~exist('no_dims', 'var') || isempty(no_dims)
    no_dims = 2;
  end
  if ~exist('perplexity', 'var') || isempty(perplexity)
    perplexity = 30;
  end
  
  % First check whether we already have an initial solution
  if numel(no_dims) > 1
    initial_solution = true;
    ydata = no_dims;
    no_dims = size(ydata, 2);
  else
    initial_solution = false;
  end
  
  % Compute joint probabilities
  D = D / max(D(:));                                                      % normalize distances
  P = d2p(D .^ 2, perplexity, 1e-5, verbose);                                      % compute affinities using fixed perplexity
  
  % Run t-SNE
  if initial_solution
    ydata = tsne_p(P, labels, ydata);
  else
    ydata = tsne_p(P, labels, no_dims);
  end
end


function [P, beta] = x2p(X, u, tol)
  %X2P Identifies appropriate sigma's to get kk NNs up to some tolerance
  %
  %   [P, beta] = x2p(xx, kk, tol)
  %
  % Identifies the required precision (= 1 / variance^2) to obtain a Gaussian
  % kernel with a certain uncertainty for every datapoint. The desired
  % uncertainty can be specified through the perplexity u (default = 15). The
  % desired perplexity is obtained up to some tolerance that can be specified
  % by tol (default = 1e-4).
  % The function returns the final Gaussian kernel in P, as well as the
  % employed precisions per instance in beta.
  %
  %
  % (C) Laurens van der Maaten, 2008
  % Maastricht University
  
  
  if ~exist('u', 'var') || isempty(u)
    u = 15;
  end
  if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
  end
  
  % Initialize some variables
  n = size(X, 1);                     % number of instances
  P = zeros(n, n);                    % empty probability matrix
  beta = ones(n, 1);                  % empty precision vector
  logU = log(u);                      % log of perplexity (= entropy)
  
  % Compute pairwise distances
  if verbose, disp('Computing pairwise distances...'); end
  sum_X = sum(X .^ 2, 2);
  D = bsxfun(@plus, sum_X, bsxfun(@plus, sum_X', -2 * (X * X')));
  
  % Run over all datapoints
  if verbose, disp('Computing P-values...'); end
  for i=1:n
    
    if verbose && ~rem(i, 500)
       disp(['Computed P-values ' num2str(i) ' of ' num2str(n) ' datapoints...']);
    end
    
    % Set minimum and maximum values for precision
    betamin = -Inf;
    betamax = Inf;
    
    % Compute the Gaussian kernel and entropy for the current precision
    Di = D(i, [1:i-1 i+1:end]);
    [H, thisP] = Hbeta(Di, beta(i));
    
    % Evaluate whether the perplexity is within tolerance
    Hdiff = H - logU;
    tries = 0;
    while abs(Hdiff) > tol && tries < 50
      
      % If not, increase or decrease precision
      if Hdiff > 0
        betamin = beta(i);
        if isinf(betamax)
          beta(i) = beta(i) * 2;
        else
          beta(i) = (beta(i) + betamax) / 2;
        end
      else
        betamax = beta(i);
        if isinf(betamin)
          beta(i) = beta(i) / 2;
        else
          beta(i) = (beta(i) + betamin) / 2;
        end
      end
      
      % Recompute the values
      [H, thisP] = Hbeta(Di, beta(i));
      Hdiff = H - logU;
      tries = tries + 1;
    end
    
    % Set the final row of P
    P(i, [1:i - 1, i + 1:end]) = thisP;
  end
  % syaoming's modification
  % disp(['Mean value of sigma: ' num2str(mean(sqrt(1 ./ beta)))]);
  % disp(['Minimum value of sigma: ' num2str(min(sqrt(1 ./ beta)))]);
  % disp(['Maximum value of sigma: ' num2str(max(sqrt(1 ./ beta)))]);
end


function [P, beta] = d2p(D, u, tol, verbose)
  %D2P Identifies appropriate sigma's to get kk NNs up to some tolerance
  %
  %   [P, beta] = d2p(D, kk, tol)
  %
  % Identifies the required precision (= 1 / variance^2) to obtain a Gaussian
  % kernel with a certain uncertainty for every datapoint. The desired
  % uncertainty can be specified through the perplexity u (default = 15). The
  % desired perplexity is obtained up to some tolerance that can be specified
  % by tol (default = 1e-4).
  % The function returns the final Gaussian kernel in P, as well as the
  % employed precisions per instance in beta.
  %
  %
  % (C) Laurens van der Maaten, 2008
  % Maastricht University
  
  
  if ~exist('u', 'var') || isempty(u)
    u = 15;
  end
  if ~exist('tol', 'var') || isempty(tol)
    tol = 1e-4;
  end
  
  % Initialize some variables
  n = size(D, 1);                     % number of instances
  P = zeros(n, n);                    % empty probability matrix
  beta = ones(n, 1);                  % empty precision vector
  logU = log(u);                      % log of perplexity (= entropy)
  
  % Run over all datapoints
  for i=1:n
    
    if verbose && ~rem(i, 500)
      % disp(['Computed P-values ' num2str(i) ' of ' num2str(n) ' datapoints...']);
    end
    
    % Set minimum and maximum values for precision
    betamin = -Inf;
    betamax = Inf;
    
    % Compute the Gaussian kernel and entropy for the current precision
    [H, thisP] = Hbeta(D(i, [1:i - 1, i + 1:end]), beta(i));
    
    % Evaluate whether the perplexity is within tolerance
    Hdiff = H - logU;
    tries = 0;
    while abs(Hdiff) > tol && tries < 50
      
      % If not, increase or decrease precision
      if Hdiff > 0
        betamin = beta(i);
        if isinf(betamax)
          beta(i) = beta(i) * 2;
        else
          beta(i) = (beta(i) + betamax) / 2;
        end
      else
        betamax = beta(i);
        if isinf(betamin)
          beta(i) = beta(i) / 2;
        else
          beta(i) = (beta(i) + betamin) / 2;
        end
      end
      
      % Recompute the values
      [H, thisP] = Hbeta(D(i, [1:i - 1, i + 1:end]), beta(i));
      Hdiff = H - logU;
      tries = tries + 1;
    end
    
    % Set the final row of P
    P(i, [1:i - 1, i + 1:end]) = thisP;
  end
  % syaoming's modification
  % disp(['Mean value of sigma: ' num2str(mean(sqrt(1 ./ beta)))]);
  % disp(['Minimum value of sigma: ' num2str(min(sqrt(1 ./ beta)))]);
  % disp(['Maximum value of sigma: ' num2str(max(sqrt(1 ./ beta)))]);
end


% Function that computes the Gaussian kernel values given a vector of
% squared Euclidean distances, and the precision of the Gaussian kernel.
% The function also computes the perplexity of the distribution.
function [H, P] = Hbeta(D, beta)
  P = exp(-D * beta);
  sumP = sum(P);
  H = log(sumP) + beta * sum(D .* P) / sumP;
  % why not: H = exp(-sum(P(P > 1e-5) .* log(P(P > 1e-5)))); ???
  P = P / sumP;
end