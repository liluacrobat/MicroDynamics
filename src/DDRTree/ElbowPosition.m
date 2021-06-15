function Para_opt = ElbowPosition(Obj, Length, parameter_range,plot_flag, EstimateSigma,num)

%% Mixture of Regression to find the two straight line
beta = RegMixEM(Length,Obj,2);
bestLength = -(beta(2,1) - beta(2,2)) / (beta(1,1) - beta(1,2));
[~,optIndex] = min(abs(Length - bestLength));
Para_opt = parameter_range(optIndex);


%% Plot

if(plot_flag)
    figure;hold on
    %plot(Length,Obj,'b*');%, 'MarkerSize', 10);

    plot([bestLength,bestLength],[0,max(Obj)*1.5],'r--','Linewidth', 3);
    if EstimateSigma==1; 
        text(Length(optIndex)+max(Length)/20,Obj(optIndex)+num,...
            ['Sigma = ' num2str(Para_opt)],'FontSize',14);
    else
        text(Length(optIndex)+max(Length)/20,Obj(optIndex)+num,...
            ['Lambda = ' num2str(Para_opt)],'FontSize',14);
    end
    xlim([0,max(Length)])
    ylim([min(Obj),max(Obj)])
    xlabel('Tree Length')
    ylabel('Mean Square Error')
    axis square;
    box on
%     boldify2
    for j = 1:2
        xplot = linspace(0,max(Length)*2,100);
        plot(xplot,xplot*beta(1,j)+beta(2,j),'b-', 'Linewidth',3);
    end
    plot(Length,Obj,'ko', 'Linewidth', 2, 'MarkerSize', 10);%, 'Color', [0.5, 0.5, 0.5])
    plot(Length(optIndex),Obj(optIndex),'ro', 'MarkerSize', 20, 'Linewidth', 3);
    hold off;
end
set(gca,'FontSize',14);
end

