function plotCellLineageLandscape(S)
% Don't use this function: 
%   some of the variable in S is discarded, thus this function is not working
%   now.

disp('plotCellLineageLandscape is invalid now due to change in S.');
return

nBranches = S.nBranches;

% dimension of subplot grids, for branch-specific plotting
n = ceil(nBranches^0.5);
m = ceil(nBranches/n);

[f,xi] = ksdensity(S.pseudotime, 0:0.01:1);
plot(xi,-log(f),'LineWidth',2)
title('Cell Lineage Landscape');
xlabel('Pseudotime');
ylabel('Energy');

for i = 1:nBranches
    branchName = ['branch' num2str(i)];
    pt = S.(branchName).pseudotime';

    subplot(m, n, i);
%     histogram(pt);
    [f,xi] = ksdensity(pt,0:0.01:1);
    plot(xi,-log(f),'LineWidth',2)
    title(sprintf('Landscape of Branch %i', i));
    xlabel('Pseudotime');
    ylabel('Energy');
end  % subplot for each branch