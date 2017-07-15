function plotNetwork( adj)
imagesc(adj,'AlphaData',0.85); colormap(flipud(gray))

NumTicks = size(adj,2)+1;
L = get(gca,'XLim');
set(gca,'XTick',linspace(L(1),L(2),NumTicks))
grid on
set(gca,'XTick',linspace(L(1)-0.5,L(2)-0.5,NumTicks))
names = 1:NumTicks-1;

stringy = {' '};
for k = 1:length(names)
    stringy = [stringy ; num2str(k)];
end


set(gca,'XTickLabel',stringy)
L = get(gca,'YLim');
set(gca,'YTick',linspace(L(1),L(2),NumTicks))
set(gca,'YTickLabel',stringy)
ylabel('node - target', 'FontSize',18);
xlabel('node - source', 'FontSize',18);
axis square


end

