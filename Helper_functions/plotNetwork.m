function plotNetwork( adj)
imagesc(adj); colormap(flipud(gray))


NumTicks = size(adj,2)+1;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
names = 1:NumTicks-1;

stringy = {' '};
for k = 1:length(names)
    stringy = [stringy ; num2str(k)];
end


set(gca,'XTickLabel',stringy)
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
set(gca,'YTickLabel',stringy)
ylabel('electrode - target', 'FontSize',14);
xlabel('electrode - source', 'FontSize',14);
axis square
grid on

end

