function [deg_ind, score] = DEG_find( expr, cutoff )
%DEGFILTERING Summary of this function goes here
%   Detailed explanation goes here

num_bins = 20;
log_mean_gene = log(mean(expr)+1);
log_dispersion_gene = log(var(expr)./mean(expr));
bins = linspace(min(log_mean_gene), max(log_mean_gene), num_bins+1);
bins(end) = bins(end) + eps;
deg_ind = [];
score = [];
for i = 1:num_bins
    index_bins = find(log_mean_gene >= bins(i) & ...
        log_mean_gene < bins(i+1));
    bin_score = zscore(log_dispersion_gene(index_bins));
    sub_index_bins = find(bin_score > cutoff );
    deg_ind = [deg_ind index_bins(sub_index_bins)];
    score = [score bin_score(sub_index_bins)];
end

end

