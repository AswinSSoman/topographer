function [toydata, initialCells] = testData(dataType,dataNumber)
n = dataNumber;
switch dataType
    case 'S'
        t = linspace(pi,0,n);
        x = cos(t);
        y = sin(2 * t);
        data = [x; y]' + randn(n, 2) / 10;
        
    case 'Y'
        n_seg = round(n/2);
        n_circle = n - n_seg;
        
        x = linspace(-1, 1, n_seg); y = exp(x);
        seg = [x; y]';
        
        t = linspace(pi*3/4, pi*9/4, n_circle);
        x = cos(t) + x(end) + 1;
        y = sin(t) + y(end);
        circle = [x; y]';
        
        data = [seg; circle] + randn(n, 2) / 5;
        
    case 'DoubleBranch'
        n_line = round(n/7);
        n_circle = round(n/3);
        n_wave = n - 2 * n_circle - n_line;
        
        x = linspace(0, 2*pi, n_wave); y = sin(x) ;
        wave = [x; y]';
        
        n_up = round(n_circle/3); n_down = n_circle - n_up;
        t = [ linspace(pi/2, pi, n_up)  linspace(pi, pi*3/2, n_down) ];
        x = pi * cos(t); y = pi * sin(t);
        left_circle = [3*pi + x; y]';
        t_circle = [4*pi + x; y - pi]';
        
        x = linspace(pi, 3*pi, n_line); y = sin(x)*0.2 + pi + 1;
        line = [x; y]';
        
        data = [wave; left_circle; t_circle; line ] + randn(n, 2) / 2;
        
    case 'Spiral'
        t = linspace(1, 3 * pi, n);
        x = cos(t); y = sin(t); z = 2 * t;
        data = [x' y' z'] + randn(n, 3)/3;
        
    case 'Wadd' % mimicing the Waddington's illustration
        n2 = round(n / 4);
        n3 = round(n / 8);
        n4 = round(n / 5.5);
        n5 = round(n / 5.5);
        n1 = n - n2 - n3 - n4 - n5;
                
        % the stem
        spacing = [5 3 3 1 2 1 2 2 2 3 4 4 3 1 2 1 1 1 1 3 4 5];
        v_level = [1 1.2 1.3 1 1 1 2.3 2.5 2.4 2.6 2.4 1 ];
        t1 = stepspace(-2, 0.5, n1, spacing);
        x1 = t1;
        y1 = exp(t1) / 3;
        b1 = [x1 + noise1(n1, [1 1]) * 0.05; % a little noise on x 
              y1 + noise1(n1, v_level) * 0.2]';
        
        % first branch, top
        spacing = [4 2 1.9 2 1.8 2 2 2 1 1 1 1 2];
        v_level = [0.3 1 1 1 1.2 1 1.3 2 2 2 3 3 3 2];
        t2 = stepspace(2, 15, n2, spacing);
        x2 = (t2 - min(t2))/5 + x1(end);
        y2 = log(t2) - log(2) + y1(end);
        b2 = [x2 + noise1(n2, [1 1]) * 0.05; 
              y2 + noise1(n2, v_level) * 0.2]';
        
        % first branch, bottom
        spacing = [4 5 6 5 5 3 2 2 2 2 3 4 5 ];
        v_level = [2 1.6 2 2 1 3 4 4 4 4 5 3];
        t3 = stepspace(-0.7 * pi, 1.5 * pi, n3, spacing);
        x3 = t3 - t3(1); x3 = x3 / 6 + x1(end);
        y3 = tanh(-t3); y3 = (y3 - y3(1)) / 6 + y1(end);
        b3 = [x3 + noise1(n3, [1 1]) * 0.05;
              y3 + noise1(n3, v_level) * 0.2]';

        % second branch, top
        spacing = exp(linspace(1, -2, 100));
        v_level = [1 2 2 2.3 5 6 5];
        t4 = stepspace(-1, 2, n4, spacing);
        x4 = (t4 - t4(1))/1.5 + x3(end);
        y4 = exp(t4); y4 = (y4 - y4(1)) / 10  + y3(end);
        b4 = [x4 + noise1(n4, [1 1]) * 0.07;
              y4 + noise1(n4, v_level) * 0.15]';
        
        % second branch, bottom
        spacing = [5 5 4 3 1 2];
        v_level = [1 2 2 2.3 4 5 4];
        t5 = stepspace(1, 5, n5, spacing);
        x5 = (t5 - t5(1))/2 + x3(end);
        y5 = sin(t5); y5 = (y5 - y5(1)) / 4 + y3(end);
        b5 = [x5 + noise1(n5, [1 1]) * 0.07;
              y5 + noise1(n5, v_level) * 0.15]';
        
        traj = [b1; b2; b3; b4; b5];
        data = traj;
        
    case 'Triple' % a triple-branching demonstration
        n2 = round(n * 0.25);
        n3 = round(n * 0.25);
        n4 = round(n * 0.25);
        n1 = n - n2 - n3 - n4;
        
        % the stem
        t1 = stepspace(-1.5, 0, n1, [6 5 3 1 1 1 3 4 6 7 5 4 5 3 2 1 1 2 4 5 4 9]);
%         x1 = t1; y1 = t1; z1 = t1;
        x1 = -2*exp(2*t1); y1 = -4*exp(t1); z1 = -2*exp(t1);
        b1 = [x1; y1; z1;]';
        
        shared_spacing = [6 4 3 3 2 2 1 2 3];
        
        % first branch, close to the x-axis
        t2 = stepspace(-2, -1, n2, shared_spacing);
        x2 = t2 - t2(1); x2 = x2 / range(x2) * 5; % 0~5
        y2 = exp(t2); y2 = y2 - y2(1);
        z2 = exp(t2); z2 = z2 - z2(1);
        b2 = [x2; y2; z2;]';
        
        % second branch, close to the y-axis
        t3 = stepspace(0, 2.5 * pi, n3, shared_spacing);
        x3 = zeros(1, n3);
        y3 = t3 - t3(1); y3 = y3 / range(y3) * 5; % 0~5
        z3 = sin(t3) + y3/3; z3 = z3 / 5;
        b3 = [x3; y3; z3]';
        
        % thirid branch, close to z-aixs
        t4 = stepspace(0, 2, n4, shared_spacing);
        x4 = zeros(1, n4) + 0.3 * t4;
        y4 = zeros(1, n4) + 0.3 * t4;
        z4 = t4; z4 = z4 / range(z4) * 5; % 0~5
        b4 = [x4; y4; z4]';
        
%         data = [b1; b2; b3; b4] + randn(n, 3) / 6;
        data = [b1 + randn(n1,3)/5; b2 + randn(n2,3)/4; b3 + randn(n3,3)/4; b4 + randn(n4,3)/5];
end

p = randperm(n);
toydata = data(p, :);
initialCells = find(p <= 3); % index of the first three cells
end

function x = stepspace(a, b, n, spacing)
assert(~any(spacing < 0));
steps = cumsum(spacing);
fullsteps = quantile(steps, n);
x = (fullsteps - min(fullsteps)) * ( (b-a) / range(steps)) + a;
end

function v = level_expand(n, level_sample, up_limit)
% expand and reample the given level_sample array 
if ~exist('up_limit', 'var'), up_limit = max(level_sample); end
assert(all(level_sample <= up_limit));

m = numel(level_sample);
pad = min(level_sample);

level_sample_padded = [pad level_sample pad]; % padding 
q_points = ( (1:n) - 0.5 ) / n * (m + 1) + 1; % resample like in image 
v = interp1(level_sample_padded, q_points) / up_limit;
end

function ns = noise1(n, noise_level)
src = min(randn(n, 1), 1.5); % preven extreme value 
ns = src' .* level_expand(n, noise_level);
end
function  ns = noise2(n, noise_level) 
ns = randn(n, 2) .* repmat(level_expand(n, noise_level)', 1, 2);
end

