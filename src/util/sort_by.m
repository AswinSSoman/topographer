function [ sorted_data ] = sort_by( data, compare_info, mode )
%SORT_BY Summary of this function goes here
%   Detailed explanation goes here

if ~exist('mode', 'var'), mode = 'ascend'; end

[~, order] = sort(compare_info, mode);
sorted_data = data(order);

end
