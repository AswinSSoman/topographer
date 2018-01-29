function [ pts, pt_indices, pt_mean_intensity, pt_std_intensity ] = ...
    wanderlust_line_graph( wanderlust, data, legend_labels, num_pts, ...
    window_width, v )
% [ pts, pt_indices, pt_mean_intensity, pt_std_intensity ] = 
%    wanderlust_line_graph( wanderlust, data, markers, num_pts, window_width )
%
% plot a line graph of marker intensities across the trajectory.
%
% (for each of num_pts equidistant points across wanderlust, pick all points in
% window_width
% diameter, use these to plot)
%  
% TODO: rewrite this travesty of a comment.

group_ids = unique(v);
if (length(group_ids) > 1) 
    legend_labels = remove_repeating_strings(legend_labels);
    legend_labels{end+1} = 'average';
end

% normalize each marker to 0-1 range using top and bottom one percentiles
% norm_bottom = prctile( data, 2 );
% norm_top = prctile( data, 98 );
% data = ( data - repmat( norm_bottom, length( data ), 1 ) ) ./ ( repmat( norm_top, length( data ), 1 ) - repmat( norm_bottom, length( data ), 1 ) );
% data = max(data, 0);
% data = min(data, 1);

%normalize the wanderlust channel?? but whyyyyhhhhaaahhhyyyy?
norm_bottom = prctile( wanderlust, 2 );
norm_top = prctile( wanderlust, 98 );
wanderlust(wanderlust<norm_bottom) = norm_bottom;
wanderlust(wanderlust>norm_top) = norm_top;
wanderlust = (wanderlust-min(wanderlust));
wanderlust = wanderlust ./ max(wanderlust);


% start from num_pts points
pts = linspace( min( wanderlust ), max( wanderlust ), num_pts );

% find all points in window_width diameter around each pts
for pt_idx = 1:num_pts
	pt = pts( pt_idx );
	pt_indices{pt_idx} = find( pt - window_width / 2 < wanderlust & wanderlust < pt + window_width / 2 );
end
 
% for each point, calculate mean of each marker around that point
pt_mean_intensity = zeros(num_pts, length(legend_labels));
pt_std_intensity  = zeros(num_pts, length(legend_labels));

if size( data, 2 )>1
    for pt_idx = 1:num_pts
        for marker_idx = 1:size( data, 2 )
            d = data( pt_indices{pt_idx}, marker_idx );
            pt_mean_intensity( pt_idx, marker_idx ) = mean( d );
            pt_std_intensity( pt_idx, marker_idx ) = std( d );
        end
    end
else % calc by grouping
    for pt_idx = 1:num_pts
        d = data( pt_indices{pt_idx});

        pt_mean_intensity( pt_idx, end ) = mean( d );
        pt_std_intensity( pt_idx, end ) = std( d );
        
        % if more than one grouping is specified, calc the mean for each
        if (length(group_ids) > 1)
            groupings = v( pt_indices{pt_idx} );

            pt_mean_intensity( pt_idx, 1:length(group_ids) ) = arrayfun(@(k) mean( d(groupings==k)), group_ids);% ./ mean( d );
            pt_std_intensity( pt_idx, 1:length(group_ids) ) = arrayfun(@(k) std( d(groupings==k)), group_ids);% ./ std( d );
        end
    end
end

% plot
colors = distinguishable_colors( size( pt_mean_intensity, 2 ) );

% shift so the generic blue color is last 
colors = [colors(2:end, :); colors(1, :)];

hold on;
cols = 1:size(pt_mean_intensity,2);
for col = cols
	plot( pts, pt_mean_intensity( :, col ), 'LineWidth', 2, 'Color', colors( col, : ));
end
hold off;


% embellish
% axis( [ 0 1 0 1 ] );
view(2);
if length(legend_labels) > 1 % syaoming's modification
  legend( legend_labels, 'Location', 'NorthEastOutside', 'Interpreter', 'None' );
end
graphLabels( '', 'Trajectory', 'Marker intensity (normalized)' );


end


function colors = distinguishable_colors(n_colors,bg,func)
% DISTINGUISHABLE_COLORS: pick colors that are maximally perceptually distinct
%
% When plotting a set of lines, you may want to distinguish them by color.
% By default, Matlab chooses a small set of colors and cycles among them,
% and so if you have more than a few lines there will be confusion about
% which line is which. To fix this problem, one would want to be able to
% pick a much larger set of distinct colors, where the number of colors
% equals or exceeds the number of lines you want to plot. Because our
% ability to distinguish among colors has limits, one should choose these
% colors to be "maximally perceptually distinguishable."
%
% This function generates a set of colors which are distinguishable
% by reference to the "Lab" color space, which more closely matches
% human color perception than RGB. Given an initial large list of possible
% colors, it iteratively chooses the entry in the list that is farthest (in
% Lab space) from all previously-chosen entries. While this "greedy"
% algorithm does not yield a global maximum, it is simple and efficient.
% Moreover, the sequence of colors is consistent no matter how many you
% request, which facilitates the users' ability to learn the color order
% and avoids major changes in the appearance of plots when adding or
% removing lines.
%
% Syntax:
%   colors = distinguishable_colors(n_colors)
% Specify the number of colors you want as a scalar, n_colors. This will
% generate an n_colors-by-3 matrix, each row representing an RGB
% color triple. If you don't precisely know how many you will need in
% advance, there is no harm (other than execution time) in specifying
% slightly more than you think you will need.
%
%   colors = distinguishable_colors(n_colors,bg)
% This syntax allows you to specify the background color, to make sure that
% your colors are also distinguishable from the background. Default value
% is white. bg may be specified as an RGB triple or as one of the standard
% "ColorSpec" strings. You can even specify multiple colors:
%     bg = {'w','k'}
% or
%     bg = [1 1 1; 0 0 0]
% will only produce colors that are distinguishable from both white and
% black.
%
%   colors = distinguishable_colors(n_colors,bg,rgb2labfunc)
% By default, distinguishable_colors uses the image processing toolbox's
% color conversion functions makecform and applycform. Alternatively, you
% can supply your own color conversion function.
%
% Example:
%   c = distinguishable_colors(25);
%   figure
%   image(reshape(c,[1 size(c)]))
%
% Example using the file exchange's 'colorspace':
%   func = @(x) colorspace('RGB->Lab',x);
%   c = distinguishable_colors(25,'w',func);

% Copyright 2010-2011 by Timothy E. Holy

  % Parse the inputs
  if (nargin < 2)
    bg = [1 1 1];  % default white background
  else
    if iscell(bg)
      % User specified a list of colors as a cell aray
      bgc = bg;
      for i = 1:length(bgc)
	bgc{i} = parsecolor(bgc{i});
      end
      bg = cat(1,bgc{:});
    else
      % User specified a numeric array of colors (n-by-3)
      bg = parsecolor(bg);
    end
  end
  
  % Generate a sizable number of RGB triples. This represents our space of
  % possible choices. By starting in RGB space, we ensure that all of the
  % colors can be generated by the monitor.
  n_grid = 30;  % number of grid divisions along each axis in RGB space
  x = linspace(0,1,n_grid);
  [R,G,B] = ndgrid(x,x,x);
  rgb = [R(:) G(:) B(:)];
  if (n_colors > size(rgb,1)/3)
    error('You can''t readily distinguish that many colors');
  end
  
  % Convert to Lab color space, which more closely represents human
  % perception
  if (nargin > 2)
    lab = func(rgb);
    bglab = func(bg);
  else
    C = makecform('srgb2lab');
    lab = applycform(rgb,C);
    bglab = applycform(bg,C);
  end

  % If the user specified multiple background colors, compute distances
  % from the candidate colors to the background colors
  mindist2 = inf(size(rgb,1),1);
  for i = 1:size(bglab,1)-1
    dX = bsxfun(@minus,lab,bglab(i,:)); % displacement all colors from bg
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
  end
  
  % Iteratively pick the color that maximizes the distance to the nearest
  % already-picked color
  colors = zeros(n_colors,3);
  lastlab = bglab(end,:);   % initialize by making the "previous" color equal to background
  for i = 1:n_colors
    dX = bsxfun(@minus,lab,lastlab); % displacement of last from all colors on list
    dist2 = sum(dX.^2,2);  % square distance
    mindist2 = min(dist2,mindist2);  % dist2 to closest previously-chosen color
    [~,index] = max(mindist2);  % find the entry farthest from all previously-chosen colors
    colors(i,:) = rgb(index,:);  % save for output
    lastlab = lab(index,:);  % prepare for next iteration
  end
    colors = min(.8, colors+.15); % make pretty
end

function c = parsecolor(s)
  if ischar(s)
    c = colorstr2rgb(s);
  elseif isnumeric(s) && size(s,2) == 3
    c = s;
  else
    error('MATLAB:InvalidColorSpec','Color specification cannot be parsed.');
  end
end

function c = colorstr2rgb(c)
  % Convert a color string to an RGB value.
  % This is cribbed from Matlab's whitebg function.
  % Why don't they make this a stand-alone function?
  rgbspec = [1 0 0;0 1 0;0 0 1;1 1 1;0 1 1;1 0 1;1 1 0;0 0 0];
  cspec = 'rgbwcmyk';
  k = find(cspec==c(1));
  if isempty(k)
    error('MATLAB:InvalidColorString','Unknown color string.');
  end
  if k~=3 || length(c)==1,
    c = rgbspec(k,:);
  elseif length(c)>2,
    if strcmpi(c(1:3),'bla')
      c = [0 0 0];
    elseif strcmpi(c(1:3),'blu')
      c = [0 0 1];
    else
      error('MATLAB:UnknownColorString', 'Unknown color string.');
    end
  end
end


function graphLabels( titleLabel, xlabelLabel, ylabelLabel, fontSize, fontWeight )
% graphLabels( titleLabel, xlabelLabel, ylabelLabel, fontSize = 14, fontWeight = 'bold' )
% sets up graph labels according to my favorite style (size 14, bold).

if( nargin < 5 )
	fontWeight = 'bold';
end

if( nargin < 4 )
	fontSize = 14;
end

if( nargin >= 1 )
	title( titleLabel, 'FontSize', fontSize, 'FontWeight', fontWeight )
end

if( nargin >= 2 )
	xlabel( xlabelLabel, 'FontSize', fontSize, 'FontWeight', fontWeight )
end

if( nargin >= 3 )
	ylabel( ylabelLabel, 'FontSize', fontSize, 'FontWeight', fontWeight )
end

end

function result = remove_repeating_strings(strings)
    if isempty(strings)
        return;
    end
    
    if numel(strings) == 1
        result = strings;
        return;
    end
    
    str = strings{1};
    for i=1:numel(str)
        TF = strncmpi(str(1:i),strings,i);
        if ~isempty(find(TF==0))
            break;
        end
    end
    
    if i>1
        result = cellfun(@(str)(str(i:end)), strings, 'UniformOutput', false);
    else
        result = strings;
    end
end
