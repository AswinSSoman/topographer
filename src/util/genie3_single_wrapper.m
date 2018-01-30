function weights = genie3_single_wrapper(expr_data, target_idx, source_idx)
% Wrap GENIE3 function to ignore the printings.

if ~exist('source_idx', 'var'), source_idx = 1:size(expr_data, 2); end;
assert(numel(target_idx) == 1);
assert(numel(source_idx) > 0);
assert(numel(source_idx) <= size(expr_data, 2));
assert(max(source_idx) <= size(expr_data, 2));

% NOTE: regulation_to_targer has zeros value for non-source_idx
cmd = 'genie3_single(expr_data, target_idx, source_idx);';
[~, regulation_to_targer] = evalc(cmd); % suppress the annoying prints
weights = regulation_to_targer(source_idx);

end

