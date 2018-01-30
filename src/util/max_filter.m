function [ ret ] = max_filter( x, window_size )
%MAX_FILTER Summary of this function goes here
%   Detailed explanation goes here

n = numel(x);
ret = zeros(n, 1);

r = floor(window_size / 2);
for i = 1:n
    i0 = max(1, i - r);
    i1 = min(n, i + r);
    ret(i) = max(x(i0:i1));
end

end

