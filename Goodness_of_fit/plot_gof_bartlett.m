function plot_gof_bartlett(result, Ij)

  n_sensors = size(result,1);    

  for nk=1:n_sensors
    figure
    %%%% Plot the data, and the estimate.
    subplot(3,1,1)
    plot(result{nk}.t,result{nk}.data, 'k')
    hold on
    plot(result{nk}.t,result{nk}.estimate, 'r')
    hold off
    axis tight
    xlabel('Time [s]')
    title(['Bartlett sensor ' num2str(nk)'])
  
    subplot(3,1,2)
    loglog(result{nk}.fj, result{nk}.f/var(result{nk}.data), 'k')
    hold on
    plot(result{nk}.fj, result{nk}.I_NX/var(result{nk}.estimate), 'r')
    plot([min(Ij(:,2)), min(Ij(:,2))], [min(result{nk}.I_NX), max(result{nk}.I_NX)], 'b:')
    plot([max(Ij(:,2)), max(Ij(:,2))], [min(result{nk}.I_NX), max(result{nk}.I_NX)], 'b:')
    hold off
    xlabel('Frequency [Hz]')
    ylabel('Power (normalized)')
    axis tight
    
    subplot(3,1,3)
    x = (0:0.1:40);
    chi2 = zeros(size(x));
    for j=1:length(x)
        chi2(j) = chi2pdf(x(j),result{nk}.k-1);
    end
    plot(x,chi2)
    hold on
    plot([result{nk}.M,result{nk}.M], [0, max(chi2)], 'r:')
    hold off
    xlabel('M values')
    ylabel('chi2')
    title(['p value is ' num2str(result{nk}.p,3)])
  end
end

