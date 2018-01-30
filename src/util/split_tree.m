function [ complete_branches ] = split_tree( ltree, subname, field_name )
%SPLIT_TREE concates a specific fiels of the tree nodes into complete branches,
% each branches may be identical in most of the heading entries. 
%
% OUTPUT: 
% complete_branches: is an (vertical) array. Each cell reprents a branch. 
% 

assert(ischar(field_name));

iter = 1;
complete_branches = {};
while isfield(ltree, subname(iter))
    subtree = ltree.(subname(iter));
    
    % recursive calls
    sub_braches_a = split_tree( subtree, subname, field_name);
    
    for b = 1:numel(sub_braches_a)
        subbranch = sub_braches_a{b};
        
        % concate a sub branch
        subbranch = [ ltree.(field_name) subbranch ]; 
        
        % add to the return 
        complete_branches = [ complete_branches ; subbranch ];
    end
    
    iter = iter + 1;
end

% if no sub-trees, just return this branch 
if iter == 1
    complete_branches = { ltree.(field_name) };
end


end

