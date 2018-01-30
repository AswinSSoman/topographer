function fig = save_plot(filename, fig, other_formats)
% save_plot configure the figure according to some difined standards, and save
% the plot as .fig file(s) in ./figures. 
%
% INPUT: 
%   filename: the name of the target file. 
%   fig: figure handle of the target plot. Default to be current figure. 
%   other_format: an array of formats, like {'png', 'svg'}. A copy of the plot
%       will be saved along side the default .fig file. 

if ~exist('fig', 'var') || isempty(fig), fig = gcf; end

% set size
set(fig, 'PaperUnits', 'centimeters');
allAxes = fig.Children;
if numel(allAxes) > 4
    set(fig, 'PaperPosition', [0 0 17.35 10.8]);
else
    set(fig, 'PaperPosition', [0 0 10 7.5]);
end

% set font
% for a = allAxes(:)' % avaliable for each axes
%     if ishandle(a)
%         set(a, 'FontUnits', 'points', 'FontSize', 10);
%         set(a, 'FontName', 'Helvetica', 'FontWeight', 'bold');
%     end
% end

% save files
savefig(fig, ['.\figures\' filename]);

if exist('other_formats', 'var')
    if ~iscell(other_formats), other_formats = {other_formats}; end
    fn = numel(other_formats);
    for i = 1:fn
        fmt = other_formats{i};
        saveas(fig, ['./figures/' filename], fmt);
    end
end

end