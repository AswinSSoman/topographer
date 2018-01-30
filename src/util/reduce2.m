function data_2d = reduce2(data, method, target_dim, force)
%REDUCE2 Carry out dimension reduction in two optional way, either by PCA or
% tSNE, specially for drawing (when dim = 2).

if ~exist('method', 'var'), method = 'pca'; end
if ~exist('dim', 'var'), target_dim = 2; end
if ~exist('force', 'var'), force = false; end

data_dim = size(data, 2);

if data_dim < target_dim
    fprintf('Error (reduce2): data dimension is too low (%i < %i).\n', ...
        data_dim, target_dim);
elseif (data_dim == target_dim) && ~force
    disp('Warnning(reduce2): target_dim == data_dim, reduction skipped.');
    data_2d = data;
else
    tic;
    switch method
        case 'pca'
            [~, data_2d, lt] = pca(data, 'NumComponents', target_dim);
            fprintf('pca time: %.2f sec\n', toc);
            expalained = sum(lt(1:target_dim)) / sum(lt) * 100;
            fprintf('\t lo-dim PCA explained: %.1f%% \n', expalained);
        case 'tsne'
            pca_dims = min(data_dim, 50); % tSNE applys PCA as preprocessing
            perplexity = 30; % default is 30
            labels = [];
            vb = false;
            data_2d = tsne(data, labels, target_dim, pca_dims, perplexity, vb);
            fprintf('tsne time: %.2f sec\n', toc);
    end % end switch
end % end if

end % end function

