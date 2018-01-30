function [ peaks, delta, gamma, pointer ] = density_peaks( dist, rho, k, HARD )
%DENSITY_PEAKS implements Rodriguez's method to find density peaks. 
%
% INPUT:
%   rho: local density of each object
%   dist: pair-wise distance (matrix) of all objects.
%   k: number of expected peaks. Actual number of results can be larger than this
%      quantaty. 
%
% OUTPUT:
%   peaks: indices of found density peaks.
%   delta: see Rodriguez' paper.
%   gamma: see Rodriguez' paper. 

if ~exist('HARD', 'var'), HARD = false; end

num = numel(rho);

% assert(issymmetric(dist));
% assert(num == size(dist, 1));
% assert(num >= k);

% calculate delta 
[~, rho_order] = sort(rho, 'descend');
if num > 0
    delta = nan(1, num);
    pointer = zeros(1, num);
else
    delta = [];
    pointer = [];
end

for i = 2:num  % ignore the density-maximum point 
    this = rho_order(i);
    delta(this) = inf;  % initializaed as inf
    for j = 1:(i-1)  % for those point with higher density than this
        that = rho_order(j);
        if dist(this, that) < delta(this), 
            delta(this) = dist(this, that);
            pointer(this) = that;
        end
    end
end

% special case for the density maximum (the highest peak)
if num>0
    [v, i] = max(delta);
    delta( rho_order(1) ) = v;
    pointer( rho_order(1) ) = i;
end

% assert(~any(isnan(delta)));
% assert(~any(pointer == 0));


% get gamma 
gamma = rho .* delta; 
[~, gamma_order] = sort(gamma, 'descend');

% find peaks
% find peaks
if num >= k
    rho_limit = min( rho( gamma_order(1:k) ) );
    delta_limit = min( delta( gamma_order(1:k) ) );
    peaks = find( (rho >= rho_limit) & (delta >= delta_limit) );
elseif num>0 && num<k;
    peaks = [ones(1,k-num)*gamma_order(1) gamma_order];
else
    peaks = [];
end
end

