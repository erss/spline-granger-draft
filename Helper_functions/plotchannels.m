
function plotchannels(varargin)

% PLOTCHANNELS plots each channels on a different line
% 
% PLOTCHANNELS([T,] CHANNELS) plots data of each channel on a separate
% ligne. Data should be NxM with N the number time step and M the number of
% channels. If the time axis T is not specified, it replaced by 1 : N. If T
% = [T1,...,TM] each channel is plotted with its own time axis.
% 
% PLOTCHANNELS(T, CHANNELS, CHANNELLABELS) plots with labels for channels
% 
% PLOTCHANNELS(T, CHANNELS, SUBPLOTHANDLE) plots in a given subplot
% 
% See also: PLOTMAPCHANNELS
% 
% Author: Louis-Emmnuel Martinet <louis.emmanuel.martinet@gmail.com>

FACTOR = 0.5;

if nargin == 3 && ishandle(varargin{3}(1))
    t = varargin{1};
    data = varargin{2};    
    ax = varargin{3};
    set(gcf, 'CurrentAxes', ax); 
else
    if nargin > 1
        t = varargin{1};
        data = varargin{2};
        xstr = 'time (s)';
    elseif nargin == 1
        data = varargin{1};
        t = 1 : size(data, 1);
        xstr = 'sample index';
    end
    xlim([t(1) t(end)]);
    xlabel(xstr);
    ylim([-0.5, size(data, 2) + 1.5]);
    %ylabel('channel index ');
    %set(gca, 'YTick', 1 : size(data, 2));
    if nargin > 2
      %  set(gca, 'YTickLabel', varargin{3});
        if isempty(varargin{3})
            ylabel('');
        end
    end
    ax = gca;
end
set(ax, 'NextPlot', 'replacechildren');

% Centering on 0 mean
nm = nanmean(data, 1);
data = bsxfun(@minus, data, nm);
% Rescaling factor
d = FACTOR * nanmean(max(data, [], 1) - min(data, [], 1));
% d = (max(data(:)) - min(data(:))) / 4; % bad for noisy channels
data = bsxfun(@plus, data / d, 1 : size(data, 2));
plot(t, data', varargin{4:end});
%%%%%%%% erss
set(gca,'YTickLabel',[]);
set(gca,'XTickLabel',[]);
set(gca,'XTick',[]);
set(gca,'YTick',[]);

end