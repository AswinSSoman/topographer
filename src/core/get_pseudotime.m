function S = get_pseudotime(S, data, ~)
%GETPSEUDOTIME Generate Pseudotime ordering for each cell, and project all
% cells onto the trajctory. 

% some alias
dist        = S.dist_mat;
traj_cells  = S.trajectory_cells;
pred        = S.trajectory_pred_link;
traj_assign = S.t_assign;
progen      = S.progen_cell;

% initialization
n = size(data, 1);
pdt = NaN(1, n);  % PDT = Pseu-Do-Time
pdt(progen) = 0;

% For trajectory cells, pseudotime is geodesic distance to the progenitor
pdt = getTimeRecuseive(progen, pdt, dist, pred);
% check the result, to debug
assert(all(~isnan(pdt(traj_cells))));
assert(all(~isinf(pdt(traj_cells))));

% For other cells, compute pseudotime by projection
for i = 1:n
    
    if ismember(i, traj_cells)
        continue
    end
    
    itraj = traj_assign(i); % which traj cell I am assigned to
    a = dist(i, itraj);
    
    if itraj == progen % at the root end
        
        [~, frontAngle] = getFrontPeak(i, itraj, pred, dist);
        pdt(i) = pdt(itraj) + a * cos(frontAngle);
        
    elseif ~ismember(itraj, pred) % at the leaf end
        
        [~, backAngle] = getBackPeak(i, itraj, pred, dist);
        pdt(i) = pdt(itraj) - a * cos(backAngle);
        
    else % intermidiate cell
        
        [~, frontAngle] = getFrontPeak(i, itraj, pred, dist);
        [~, backAngle] = getBackPeak(i, itraj, pred, dist);
        if min(frontAngle, backAngle) >= pi/2 % out of range
            pdt(i) = pdt(itraj);
        elseif frontAngle < backAngle
            pdt(i) = pdt(itraj) + a * cos(frontAngle);
        else
            pdt(i) = pdt(itraj) - a * cos(backAngle);
        end
        
    end % end if
    
end % end for

% check the results
assert(all(~isnan(pdt)));
assert(all(~isinf(pdt)));

% normalize and return 
S.pseudotime = (pdt-min(pdt))/(max(pdt)-min(pdt));

end


function pdt = getTimeRecuseive(this, pdt, D, pred)
nextCells = find(pred == this);
assert(~isnan(pdt(this)));
pdt(nextCells) = pdt(this) + D(this, nextCells);
for i = 1:length(nextCells)
    next = nextCells(i);
    pdt = getTimeRecuseive(next, pdt, D, pred);
end
end


function [frontPeak, frontAngle] = getFrontPeak(i, ipeak, predIdx, D)
candidate = find(predIdx == ipeak);
angles = arrayfun(@(f) acosByDist(D(i, ipeak), D(ipeak, f), D(i, f)), candidate);

[frontAngle, idx] = min(angles);
frontPeak = candidate(idx);
end


function [backPeak, backAngle] = getBackPeak(i, ipeak, predIdx, D)
backPeak = predIdx(ipeak);
backAngle = acosByDist(D(i, ipeak), D(ipeak, backPeak), D(i, backPeak));
end


function pst = projectAverage(angle, a, base, x1, x2)
p1 = a * cos(angle);
p2 = base - p1;
weight = [p2 p1];

pst = sum(weight .* [x1 x2]) / sum(weight); % weighted average by projection
end


function ang = acosByDist(a, b, c)
ang = acos( (a^2 + b^2 - c^2) / (2 * a * b) );
ang = real(ang); % some time be complex, for non-euclidean input
end