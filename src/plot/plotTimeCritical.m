function [ dState ] = plotTimeCritical( time, d_mat, k )
%PLOT Summary of this function goes here
%   Detailed explanation goes here

[time_sorted, t_order] = sort(time, 'ascend');

dTime = time_sorted(2:end) - time_sorted(1:end-1);
dState = arrayfun(@(i) d_mat(t_order(i), t_order(i+1)), 1:numel(time)-1);
ddState = abs( dState(2:end) - dState(1:end-1) );

dTime_smooth = max_filter(dTime, k)';
dState_smooth = max_filter(dState, k)';
ddState_smooth = max_filter(ddState, k)';

% Plot things
n_plot = 5;
subplot(n_plot, 1, 1);
plot(time_sorted); title('Just Time');
subplot(n_plot, 1, 2);
plot(time_sorted(1:end-1), [dState; dState_smooth]); title('dStage - Time');
subplot(n_plot, 1, 3);
plot(time_sorted(1:end-1), log(1 + dState ./ dTime_smooth)); title('lg( dStage/dTime ) - Time');
subplot(n_plot, 1, 4);
plot(time_sorted(1:end-2), [ddState; ddState_smooth]); title('ddStage - Time');
subplot(n_plot, 1, 5);
plot(time_sorted(1:end-2), ddState ./ dTime_smooth(1:end-1)); title('ddStage/dTime - Time');

end

