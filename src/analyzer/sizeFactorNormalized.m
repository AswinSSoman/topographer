function normalized_count_data = sizeFactorNormalized(count_data)
% sizeFactorNormalized  Normalization of the raw read counts with the same 
% method used byDESeq2.
% count_data is an expression matrix. Rows of count_data correspond to cells
% and columns correspond to genes.

[num_cell, num_gene] = size(count_data);
nonzero_gene = (sum(log(count_data),1)~=-Inf);
normalized_count_data = count_data(:,nonzero_gene);
normalized_count_data = normalized_count_data./repmat(geomean(normalized_count_data,1), num_cell, 1);
size_factor = median(normalized_count_data,2);
normalized_count_data = count_data./repmat(size_factor, 1, num_gene);

end