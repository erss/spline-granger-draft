function plotGraph( A )
% Plots graph for given adjacency matrix A.

if issymmetric(A)
    G = graph(A);
else
    G = digraph(A);
end

    h = plot(G,'MarkerSize',8,'LineWidth',1.5,'Layout','circle');
    h.NodeColor = [1 0 1];
  %  h.NodeLabel = {};
    axis square
    axis tight
    set(gca,'YTickLabel',[]);
    set(gca,'XTickLabel',[]);
  %  set(gca,'visible','off')
     set(gca,'XtickLabel',[],'YtickLabel',[]);
   % highlight(h,'NodeColor','k')


end

