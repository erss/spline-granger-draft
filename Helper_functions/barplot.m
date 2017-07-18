function barplot(labels,varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
n = length(varargin);

for k = 1:n
    result = varargin{k};
    mn  = mean(result);
    sem = std(result)/sqrt(length(result));
    bar(k,mn,'FaceAlpha',0.2,'EdgeColor','None')
    hold on
    plot([k k], [mn-2*sem, mn+2*sem],'k','LineWidth',2);  %This will add a vertical line.
    
    colormap(gray)
    set(gca, 'XTick', 1:n, 'XTickLabel', labels,'FontSize',18);
    xlim=get(gca,'xlim');
end

end

