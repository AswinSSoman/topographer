function plotTrajectory(S, X2d)
%PLOTBRANCH Draws the brancching trajectory as arrows, with density plot as 
% backgroud.

draw_density(X2d, S.rho);
draw_branching(X2d, S.progen_cell, S.trajectory_adj);

drawnow;

end

