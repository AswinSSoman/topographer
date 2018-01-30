function [content_array, links] = flatten_tree(ltree, subname, content_names)
%FLATTENTREE turns the result of 'linkedListTree.m' into more convinient
% formats. Use this as bulding block for grouping function.
%
% OUTPUT:
% content_array: struct containg contents of tree nodes. All nodes of the ltree 
%     are flatten as an array 'content_array'. Node contents are put in
%     different fields of the struct.
% links: un-weighted links between entries in 'content_array', representing the
%     original connection in ltree.
% 


% contents in this node
for i = 1:numel(content_names)
    iName = content_names{i};
    content_array.(iName) = { ltree.(iName) };
end
k = numel(content_array.(iName));
links.from = 1:k-1;
links.to   = 2:k;

iter = 1;
while isfield(ltree, subname(iter))
    subtree = ltree.(subname(iter));
    
    % recursive call
    [subttreeArray, stLinks] = flatten_tree(subtree, subname, content_names);

    % what to do is subtree is actually empty...
    if isempty(subttreeArray), warning('stArray empty'); end
    
    % add links
    offset = numel(content_array.(iName)); % any name is fine 
    links.from = [links.from k];
    links.to   = [links.to   offset+1];
    links.from = [links.from offset+stLinks.from];
    links.to   = [links.to   offset+stLinks.to];
    
    % append contents
    for i = 1:numel(content_names)
        iName = content_names{i};
        content_array.(iName) = [content_array.(iName) subttreeArray.(iName)];
    end
    
    iter = iter + 1;
end

end