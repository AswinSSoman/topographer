function [groupArray, groupAdjMat] = groupByTime( ltree, pdt, n, mode )
%GROUPBYTIME Group the cells and flatten the result in a more easy-to-use
%format.

% devide by time
divTree = divide_recursive(ltree, pdt, n, mode);

% flatten the result in to groups
subname = @(i) ['subtree', int2str(i)];
cname = {'cellgroups' };   % this name is used in groupByTimeRec() 

[contents, links] = flattenTree_special(divTree, subname, cname);
groupArray = contents.cellgroups;

% convert edge-list to adjacent matrix 
numSteps = numel(groupArray);
m = length(links.from);
assert(m == length(links.to));
groupAdjMat = sparse(links.from, links.to, ones(m, 1), numSteps, numSteps);
end


function divTree = divide_recursive(llt, pdt, n, mode)
% divideByTime divided the trees into a more fine-grain tree, according
% to pseudotime provided. 
%
% INPUT:
% segmentTree: a special structure returned by 'segmentTreeFromS.m'. This
%     structre must contains cell indices for each node.
% pdt: pseudotime for cells.
% n: divided the tree on n timepoints.
% mode: choose how the strategy of dividing. In 'fuzz' mode, dividing occurs on
%     each segment seperately, resulting approximate n groups for the longest 
%     branch; in 'clear' mode, dividing may cross a branch point, but result on
%     exactly n groups for the longest branch. Also the timpoints in consistent 
%     among branches
%
% OUTPUT:
% divTree: a finer-grain tree structure. Accepted by 'flatenTree.m' to provide
%     more convinient structures.

if ~exist('mode', 'var'), mode = 'fuzz'; end

member = llt.members;
timeStart = min(pdt(member));
timeEnd = max(pdt(member));

%% Divide the tree segments 

% useful functions
getMemberCA = @(lt, rt) { member( (lt<pdt(member)) & (pdt(member)<=rt) ) };
subtreeName = @(i) ['subtree' int2str(i)];

% allow two modes
switch mode
    
    case 'fuzz' % dividing points may not be what you expected
        
        % divede the current segement into groups
        nGroup = ceil(range(pdt(member)) / range(pdt) * n);
        tpoints = linspace(timeStart, timeEnd, nGroup + 1);
        cellgroups = arrayfun(getMemberCA, tpoints(1:end-1), tpoints(2:end));
        
        divTree.cellgroups = cellgroups;
        divTree.timepoints = (tpoints(1:end-1) + tpoints(2:end)) / 2;
        
        % iterate on all subtrees
        iter = 1;
        while isfield(llt, subtreeName(iter))
            subname = subtreeName(iter);
            sst = llt.(subname);
            divTree.(subname) = divide_recursive(sst, pdt, n-nGroup, 'fuzz');
            iter = iter + 1;
        end
        
    case 'clear' % dividing is strictly on n timepoints, evenly.   
        
        width = range(pdt) / (2*(n-1));
        fullstep = 2 * width;

        % divide the current segments
        t = min(pdt); % initial
        while t + width < timeStart,
            t = t + fullstep;
        end
        cellgroups = {}; % initial
        timepoints = []; % initial
        while t - width < timeEnd
            cellgroups = [cellgroups getMemberCA(t-width, t+width)];
            timepoints = [timepoints t];
            t = t + fullstep;
        end
        
        if isempty(cellgroups), warning('empty cellgroups'); end
                
        % iterate on all subtrees
        divTree = struct;
        iter = 1;
        while isfield(llt, subtreeName(iter))
            subname = subtreeName(iter);
            subDT = divide_recursive(llt.(subname), pdt, n, 'clear');
            % move the first group, if necessary
            if any(pdt(subDT.cellgroups{1}) <= t + width)
                cellgroups{end} = union(cellgroups{end}, subDT.cellgroups{1});
                subDT.cellgroups(1) = [];
                subDT.timepoints(1) = [];
            end
            if isempty(subDT.cellgroups)
                warning('time steps maybe too short.');
            end
            divTree.(subname) = subDT;
            iter = iter + 1;
        end
        
        divTree.cellgroups = cellgroups;
        divTree.timepoints = timepoints;
        
    otherwise
        
        error('unsupported mode');
        
end % end switch 

end % end function 


% this function is special coded for the main function; should be fixed someday, 
% along with the whole deisign of groupByTime 
function [groupContents, groupLinks] = flattenTree_special(tree, subname, content_names)
%FLATTENTREE turns the result of 'linkedListTree.m' into more convinient
% formats. Use this as bulding block for grouping function.
%
% OUTPUT:
% groupArray: array of cell indecise.
% timepoints: pseudotime assigned to each group.
% groupLinks: un-weighted links between each pair of groups.
% groupAdjMat: same information as groupLink, but as an adjacent matrix.

[groupContents, groupLinks] = flatRecursive(tree, subname, content_names);

end


function [contentArray, links] = flatRecursive(tree, subname, content_names)

% contents in this node
for i = 1:numel(content_names)
    iName = content_names{i};
    contentArray.(iName) = tree.(iName); % here is the speciality !!! 
end
k = numel(contentArray.(iName));
links.from = 1:k-1;
links.to   = 2:k;

iter = 1;
while isfield(tree, subname(iter))
    subtree = tree.(subname(iter));
    [subttreeArray, stLinks] = flatRecursive(subtree, subname, content_names);

    % what to do is subtree is actually empty...
    if isempty(subttreeArray), warning('stArray empty'); end
    
    % add links
    offset = numel(contentArray.(iName)); % any name is find 
    links.from = [links.from k];
    links.to   = [links.to   offset+1];
    links.from = [links.from offset+stLinks.from];
    links.to   = [links.to   offset+stLinks.to];
    
    % append contents
    for i = 1:numel(content_names)
        iName = content_names{i};
        contentArray.(iName) = [contentArray.(iName) subttreeArray.(iName)];
    end
    
    iter = iter + 1;
end

end