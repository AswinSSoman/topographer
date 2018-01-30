function plotExprLandscape(S)
% Don't use this function: 
%   some of the variable in S is discarded, thus this function is not working
%   now.
%
% future advise: re-write the function for preparing results (generateResult). 

% dimension of subplot grids, for gene-specific plotting
nGene = length(S.Opts.geneIDforPdf);
if nGene < 1, return; end
n = ceil(nGene^0.5);
m = ceil(nGene/n);


timePoints = S.exprPdf.time;
for pdfIdx = 1:nGene
    geneIdx = S.Opts.geneIDforPdf(pdfIdx);
    subplot(m, n, pdfIdx);
    
    exprPoints = S.exprPdf.smoothPdf{pdfIdx,1};
    exprPdf = S.exprPdf.smoothPdf{pdfIdx,2};
    colormap jet
    surf(timePoints, exprPoints, exprPdf', ...
        'FaceColor', 'interp', 'EdgeColor', [.3 .3 .3]);
        
    % add title
    if isfield(S, 'geneID');
        title(S.geneID{geneIdx});
    else
        title(sprintf('Gene Index = %i', geneIdx));
    end
    
    xlabel('Pseudotime');
    ylabel('Expression');
    zlabel('Probability Density');
    if exprPoints(end) > exprPoints(1), % prevent error
        ylim([exprPoints(1) exprPoints(end)])
    end
    pbaspect([10, 8, 1]);
end